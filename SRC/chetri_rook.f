      SUBROUTINE CHETRI_ROOK( UPLO, N, A, LDA, IPIV, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      COMPLEX            CONE, CZERO
      const              ONE = 1.0E+0, CONE = ( 1.0E+0, 0.0E+0 ), CZERO = ( 0.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, K, KP, KSTEP;
      REAL               AK, AKP1, D, T
      COMPLEX            AKKP1, TEMP
      // ..
      // .. External Functions ..
      bool               LSAME;
      COMPLEX            CDOTC
      // EXTERNAL LSAME, CDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CHEMV, CSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, MAX, REAL
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHETRI_ROOK', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Check that the diagonal matrix D is nonsingular.

      IF( UPPER ) THEN

         // Upper triangular storage: examine D from bottom to top

         DO 10 INFO = N, 1, -1
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.CZERO ) RETURN
   10    CONTINUE
      } else {

         // Lower triangular storage: examine D from top to bottom.

         DO 20 INFO = 1, N
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.CZERO ) RETURN
   20    CONTINUE
      END IF
      INFO = 0

      IF( UPPER ) THEN

         // Compute inv(A) from the factorization A = U*D*U**H.

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = 1
   30    CONTINUE

         // If K > N, exit from loop.

         IF( K.GT.N ) GO TO 70

         IF( IPIV( K ).GT.0 ) THEN

            // 1 x 1 diagonal block

            // Invert the diagonal block.

            A( K, K ) = ONE / REAL( A( K, K ) )

            // Compute column K of the inverse.

            IF( K.GT.1 ) THEN
               CALL CCOPY( K-1, A( 1, K ), 1, WORK, 1 )
               CALL CHEMV( UPLO, K-1, -CONE, A, LDA, WORK, 1, CZERO, A( 1, K ), 1 )                A( K, K ) = A( K, K ) - REAL( CDOTC( K-1, WORK, 1, A( 1, K ), 1 ) )
            END IF
            KSTEP = 1
         } else {

            // 2 x 2 diagonal block

            // Invert the diagonal block.

            T = ABS( A( K, K+1 ) )
            AK = REAL( A( K, K ) ) / T
            AKP1 = REAL( A( K+1, K+1 ) ) / T
            AKKP1 = A( K, K+1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K, K ) = AKP1 / D
            A( K+1, K+1 ) = AK / D
            A( K, K+1 ) = -AKKP1 / D

            // Compute columns K and K+1 of the inverse.

            IF( K.GT.1 ) THEN
               CALL CCOPY( K-1, A( 1, K ), 1, WORK, 1 )
               CALL CHEMV( UPLO, K-1, -CONE, A, LDA, WORK, 1, CZERO, A( 1, K ), 1 )                A( K, K ) = A( K, K ) - REAL( CDOTC( K-1, WORK, 1, A( 1, K ), 1 ) )                A( K, K+1 ) = A( K, K+1 ) - CDOTC( K-1, A( 1, K ), 1, A( 1, K+1 ), 1 )
               CALL CCOPY( K-1, A( 1, K+1 ), 1, WORK, 1 )
               CALL CHEMV( UPLO, K-1, -CONE, A, LDA, WORK, 1, CZERO, A( 1, K+1 ), 1 )                A( K+1, K+1 ) = A( K+1, K+1 ) - REAL( CDOTC( K-1, WORK, 1, A( 1, K+1 ), 1 ) )
            END IF
            KSTEP = 2
         END IF

         IF( KSTEP.EQ.1 ) THEN

            // Interchange rows and columns K and IPIV(K) in the leading
            // submatrix A(1:k,1:k)

            KP = IPIV( K )
            IF( KP.NE.K ) THEN

               IF( KP.GT.1 ) CALL CSWAP( KP-1, A( 1, K ), 1, A( 1, KP ), 1 )

               DO 40 J = KP + 1, K - 1
                  TEMP = CONJG( A( J, K ) )
                  A( J, K ) = CONJG( A( KP, J ) )
                  A( KP, J ) = TEMP
   40          CONTINUE

               A( KP, K ) = CONJG( A( KP, K ) )

               TEMP = A( K, K )
               A( K, K ) = A( KP, KP )
               A( KP, KP ) = TEMP
            END IF
         } else {

            // Interchange rows and columns K and K+1 with -IPIV(K) and
            // -IPIV(K+1) in the leading submatrix A(k+1:n,k+1:n)

            // (1) Interchange rows and columns K and -IPIV(K)

            KP = -IPIV( K )
            IF( KP.NE.K ) THEN

               IF( KP.GT.1 ) CALL CSWAP( KP-1, A( 1, K ), 1, A( 1, KP ), 1 )

               DO 50 J = KP + 1, K - 1
                  TEMP = CONJG( A( J, K ) )
                  A( J, K ) = CONJG( A( KP, J ) )
                  A( KP, J ) = TEMP
   50          CONTINUE

               A( KP, K ) = CONJG( A( KP, K ) )

               TEMP = A( K, K )
               A( K, K ) = A( KP, KP )
               A( KP, KP ) = TEMP

               TEMP = A( K, K+1 )
               A( K, K+1 ) = A( KP, K+1 )
               A( KP, K+1 ) = TEMP
            END IF

            // (2) Interchange rows and columns K+1 and -IPIV(K+1)

            K = K + 1
            KP = -IPIV( K )
            IF( KP.NE.K ) THEN

               IF( KP.GT.1 ) CALL CSWAP( KP-1, A( 1, K ), 1, A( 1, KP ), 1 )

               DO 60 J = KP + 1, K - 1
                  TEMP = CONJG( A( J, K ) )
                  A( J, K ) = CONJG( A( KP, J ) )
                  A( KP, J ) = TEMP
   60          CONTINUE

               A( KP, K ) = CONJG( A( KP, K ) )

               TEMP = A( K, K )
               A( K, K ) = A( KP, KP )
               A( KP, KP ) = TEMP
            END IF
         END IF

         K = K + 1
         GO TO 30
   70    CONTINUE

      } else {

         // Compute inv(A) from the factorization A = L*D*L**H.

         // K is the main loop index, decreasing from N to 1 in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = N
   80    CONTINUE

         // If K < 1, exit from loop.

         IF( K.LT.1 ) GO TO 120

         IF( IPIV( K ).GT.0 ) THEN

            // 1 x 1 diagonal block

            // Invert the diagonal block.

            A( K, K ) = ONE / REAL( A( K, K ) )

            // Compute column K of the inverse.

            IF( K.LT.N ) THEN
               CALL CCOPY( N-K, A( K+1, K ), 1, WORK, 1 )
               CALL CHEMV( UPLO, N-K, -CONE, A( K+1, K+1 ), LDA, WORK, 1, CZERO, A( K+1, K ), 1 )                A( K, K ) = A( K, K ) - REAL( CDOTC( N-K, WORK, 1, A( K+1, K ), 1 ) )
            END IF
            KSTEP = 1
         } else {

            // 2 x 2 diagonal block

            // Invert the diagonal block.

            T = ABS( A( K, K-1 ) )
            AK = REAL( A( K-1, K-1 ) ) / T
            AKP1 = REAL( A( K, K ) ) / T
            AKKP1 = A( K, K-1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K-1, K-1 ) = AKP1 / D
            A( K, K ) = AK / D
            A( K, K-1 ) = -AKKP1 / D

            // Compute columns K-1 and K of the inverse.

            IF( K.LT.N ) THEN
               CALL CCOPY( N-K, A( K+1, K ), 1, WORK, 1 )
               CALL CHEMV( UPLO, N-K, -CONE, A( K+1, K+1 ), LDA, WORK, 1, CZERO, A( K+1, K ), 1 )                A( K, K ) = A( K, K ) - REAL( CDOTC( N-K, WORK, 1, A( K+1, K ), 1 ) )                A( K, K-1 ) = A( K, K-1 ) - CDOTC( N-K, A( K+1, K ), 1, A( K+1, K-1 ), 1 )
               CALL CCOPY( N-K, A( K+1, K-1 ), 1, WORK, 1 )
               CALL CHEMV( UPLO, N-K, -CONE, A( K+1, K+1 ), LDA, WORK, 1, CZERO, A( K+1, K-1 ), 1 )                A( K-1, K-1 ) = A( K-1, K-1 ) - REAL( CDOTC( N-K, WORK, 1, A( K+1, K-1 ), 1 ) )
            END IF
            KSTEP = 2
         END IF

         IF( KSTEP.EQ.1 ) THEN

            // Interchange rows and columns K and IPIV(K) in the trailing
            // submatrix A(k:n,k:n)

            KP = IPIV( K )
            IF( KP.NE.K ) THEN

               IF( KP.LT.N ) CALL CSWAP( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 )

               DO 90 J = K + 1, KP - 1
                  TEMP = CONJG( A( J, K ) )
                  A( J, K ) = CONJG( A( KP, J ) )
                  A( KP, J ) = TEMP
   90          CONTINUE

               A( KP, K ) = CONJG( A( KP, K ) )

               TEMP = A( K, K )
               A( K, K ) = A( KP, KP )
               A( KP, KP ) = TEMP
            END IF
         } else {

            // Interchange rows and columns K and K-1 with -IPIV(K) and
            // -IPIV(K-1) in the trailing submatrix A(k-1:n,k-1:n)

            // (1) Interchange rows and columns K and -IPIV(K)

            KP = -IPIV( K )
            IF( KP.NE.K ) THEN

               IF( KP.LT.N ) CALL CSWAP( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 )

               DO 100 J = K + 1, KP - 1
                  TEMP = CONJG( A( J, K ) )
                  A( J, K ) = CONJG( A( KP, J ) )
                  A( KP, J ) = TEMP
  100         CONTINUE

               A( KP, K ) = CONJG( A( KP, K ) )

               TEMP = A( K, K )
               A( K, K ) = A( KP, KP )
               A( KP, KP ) = TEMP

               TEMP = A( K, K-1 )
               A( K, K-1 ) = A( KP, K-1 )
               A( KP, K-1 ) = TEMP
            END IF

            // (2) Interchange rows and columns K-1 and -IPIV(K-1)

            K = K - 1
            KP = -IPIV( K )
            IF( KP.NE.K ) THEN

               IF( KP.LT.N ) CALL CSWAP( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 )

               DO 110 J = K + 1, KP - 1
                  TEMP = CONJG( A( J, K ) )
                  A( J, K ) = CONJG( A( KP, J ) )
                  A( KP, J ) = TEMP
  110         CONTINUE

               A( KP, K ) = CONJG( A( KP, K ) )

               TEMP = A( K, K )
               A( K, K ) = A( KP, KP )
               A( KP, KP ) = TEMP
            END IF
         END IF

         K = K - 1
         GO TO 80
  120    CONTINUE
      END IF

      RETURN

      // End of CHETRI_ROOK

      }
