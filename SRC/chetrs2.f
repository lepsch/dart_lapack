      SUBROUTINE CHETRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = (1.0E+0,0.0E+0) )
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, IINFO, J, K, KP;
      REAL               S
      COMPLEX            AK, AKM1, AKM1K, BK, BKM1, DENOM
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CSSCAL, CSYCONV, CSWAP, CTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, MAX, REAL
      // ..
      // .. Executable Statements ..

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHETRS2', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN

      // Convert A

      CALL CSYCONV( UPLO, 'C', N, A, LDA, IPIV, WORK, IINFO )

      IF( UPPER ) THEN

         // Solve A*X = B, where A = U*D*U**H.

        // P**T * B
        K=N
        DO WHILE ( K .GE. 1 )
         IF( IPIV( K ).GT.0 ) THEN
            // 1 x 1 diagonal block
            // Interchange rows K and IPIV(K).
            KP = IPIV( K )
            IF( KP.NE.K ) CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K=K-1
         ELSE
            // 2 x 2 diagonal block
            // Interchange rows K-1 and -IPIV(K).
            KP = -IPIV( K )
            IF( KP.EQ.-IPIV( K-1 ) ) CALL CSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )
            K=K-2
         END IF
        END DO

*  Compute (U \P**T * B) -> B    [ (U \P**T * B) ]

        CALL CTRSM('L','U','N','U',N,NRHS,ONE,A,LDA,B,LDB)

*  Compute D \ B -> B   [ D \ (U \P**T * B) ]

         I=N
         DO WHILE ( I .GE. 1 )
            IF( IPIV(I) .GT. 0 ) THEN
              S = REAL( ONE ) / REAL( A( I, I ) )
              CALL CSSCAL( NRHS, S, B( I, 1 ), LDB )
            ELSEIF ( I .GT. 1) THEN
               IF ( IPIV(I-1) .EQ. IPIV(I) ) THEN
                  AKM1K = WORK(I)
                  AKM1 = A( I-1, I-1 ) / AKM1K
                  AK = A( I, I ) / CONJG( AKM1K )
                  DENOM = AKM1*AK - ONE
                  DO 15 J = 1, NRHS
                     BKM1 = B( I-1, J ) / AKM1K
                     BK = B( I, J ) / CONJG( AKM1K )
                     B( I-1, J ) = ( AK*BKM1-BK ) / DENOM
                     B( I, J ) = ( AKM1*BK-BKM1 ) / DENOM
 15              CONTINUE
               I = I - 1
               ENDIF
            ENDIF
            I = I - 1
         END DO

       // Compute (U**H \ B) -> B   [ U**H \ (D \ (U \P**T * B) ) ]

         CALL CTRSM('L','U','C','U',N,NRHS,ONE,A,LDA,B,LDB)

        // P * B  [ P * (U**H \ (D \ (U \P**T * B) )) ]

        K=1
        DO WHILE ( K .LE. N )
         IF( IPIV( K ).GT.0 ) THEN
            // 1 x 1 diagonal block
            // Interchange rows K and IPIV(K).
            KP = IPIV( K )
            IF( KP.NE.K ) CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K=K+1
         ELSE
            // 2 x 2 diagonal block
            // Interchange rows K-1 and -IPIV(K).
            KP = -IPIV( K )
            IF( K .LT. N .AND. KP.EQ.-IPIV( K+1 ) ) CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K=K+2
         ENDIF
        END DO

      ELSE

         // Solve A*X = B, where A = L*D*L**H.

        // P**T * B
        K=1
        DO WHILE ( K .LE. N )
         IF( IPIV( K ).GT.0 ) THEN
            // 1 x 1 diagonal block
            // Interchange rows K and IPIV(K).
            KP = IPIV( K )
            IF( KP.NE.K ) CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K=K+1
         ELSE
            // 2 x 2 diagonal block
            // Interchange rows K and -IPIV(K+1).
            KP = -IPIV( K+1 )
            IF( KP.EQ.-IPIV( K ) ) CALL CSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )
            K=K+2
         ENDIF
        END DO

*  Compute (L \P**T * B) -> B    [ (L \P**T * B) ]

        CALL CTRSM('L','L','N','U',N,NRHS,ONE,A,LDA,B,LDB)

*  Compute D \ B -> B   [ D \ (L \P**T * B) ]

         I=1
         DO WHILE ( I .LE. N )
            IF( IPIV(I) .GT. 0 ) THEN
              S = REAL( ONE ) / REAL( A( I, I ) )
              CALL CSSCAL( NRHS, S, B( I, 1 ), LDB )
            ELSE
                  AKM1K = WORK(I)
                  AKM1 = A( I, I ) / CONJG( AKM1K )
                  AK = A( I+1, I+1 ) / AKM1K
                  DENOM = AKM1*AK - ONE
                  DO 25 J = 1, NRHS
                     BKM1 = B( I, J ) / CONJG( AKM1K )
                     BK = B( I+1, J ) / AKM1K
                     B( I, J ) = ( AK*BKM1-BK ) / DENOM
                     B( I+1, J ) = ( AKM1*BK-BKM1 ) / DENOM
 25              CONTINUE
                  I = I + 1
            ENDIF
            I = I + 1
         END DO

*  Compute (L**H \ B) -> B   [ L**H \ (D \ (L \P**T * B) ) ]

        CALL CTRSM('L','L','C','U',N,NRHS,ONE,A,LDA,B,LDB)

        // P * B  [ P * (L**H \ (D \ (L \P**T * B) )) ]

        K=N
        DO WHILE ( K .GE. 1 )
         IF( IPIV( K ).GT.0 ) THEN
            // 1 x 1 diagonal block
            // Interchange rows K and IPIV(K).
            KP = IPIV( K )
            IF( KP.NE.K ) CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K=K-1
         ELSE
            // 2 x 2 diagonal block
            // Interchange rows K-1 and -IPIV(K).
            KP = -IPIV( K )
            IF( K.GT.1 .AND. KP.EQ.-IPIV( K-1 ) ) CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K=K-2
         ENDIF
        END DO

      END IF

      // Revert A

      CALL CSYCONV( UPLO, 'R', N, A, LDA, IPIV, WORK, IINFO )

      RETURN

      // End of CHETRS2

      }
