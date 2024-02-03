      SUBROUTINE ZHETRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16       A( LDA, * ), B( LDB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE
      const              ONE = (1.0D+0,0.0D+0) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, IINFO, J, K, KP;
      double             S;
      COMPLEX*16         AK, AKM1, AKM1K, BK, BKM1, DENOM
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZDSCAL, ZSYCONV, ZSWAP, ZTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCONJG, MAX
      // ..
      // .. Executable Statements ..

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( NRHS.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -8
      }
      if ( INFO != 0 ) {
         xerbla('ZHETRS2', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) RETURN;

      // Convert A

      zsyconv(UPLO, 'C', N, A, LDA, IPIV, WORK, IINFO );

      if ( UPPER ) {

         // Solve A*X = B, where A = U*D*U**H.

        // P**T * B
        K=N
        DO WHILE ( K .GE. 1 )
         if ( IPIV( K ).GT.0 ) {
            // 1 x 1 diagonal block
            // Interchange rows K and IPIV(K).
            KP = IPIV( K )
            if (KP != K) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            K=K-1
         } else {
            // 2 x 2 diagonal block
            // Interchange rows K-1 and -IPIV(K).
            KP = -IPIV( K )
            IF( KP == -IPIV( K-1 ) ) CALL ZSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )
            K=K-2
         }
        }

*  Compute (U \P**T * B) -> B    [ (U \P**T * B) ]

        ztrsm('L','U','N','U',N,NRHS,ONE,A,LDA,B,LDB);

*  Compute D \ B -> B   [ D \ (U \P**T * B) ]

         I=N
         DO WHILE ( I .GE. 1 )
            if ( IPIV(I) .GT. 0 ) {
              S = DBLE( ONE ) / DBLE( A( I, I ) )
              zdscal(NRHS, S, B( I, 1 ), LDB );
            } else if ( I .GT. 1) {
               if ( IPIV(I-1) == IPIV(I) ) {
                  AKM1K = WORK(I)
                  AKM1 = A( I-1, I-1 ) / AKM1K
                  AK = A( I, I ) / DCONJG( AKM1K )
                  DENOM = AKM1*AK - ONE
                  for (J = 1; J <= NRHS; J++) { // 15
                     BKM1 = B( I-1, J ) / AKM1K
                     BK = B( I, J ) / DCONJG( AKM1K )
                     B( I-1, J ) = ( AK*BKM1-BK ) / DENOM
                     B( I, J ) = ( AKM1*BK-BKM1 ) / DENOM
                 } // 15
               I = I - 1
               }
            }
            I = I - 1
         }

       // Compute (U**H \ B) -> B   [ U**H \ (D \ (U \P**T * B) ) ]

         ztrsm('L','U','C','U',N,NRHS,ONE,A,LDA,B,LDB);

        // P * B  [ P * (U**H \ (D \ (U \P**T * B) )) ]

        K=1
        DO WHILE ( K .LE. N )
         if ( IPIV( K ).GT.0 ) {
            // 1 x 1 diagonal block
            // Interchange rows K and IPIV(K).
            KP = IPIV( K )
            if (KP != K) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            K=K+1
         } else {
            // 2 x 2 diagonal block
            // Interchange rows K-1 and -IPIV(K).
            KP = -IPIV( K )
            IF( K .LT. N && KP == -IPIV( K+1 ) ) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K=K+2
         }
        }

      } else {

         // Solve A*X = B, where A = L*D*L**H.

        // P**T * B
        K=1
        DO WHILE ( K .LE. N )
         if ( IPIV( K ).GT.0 ) {
            // 1 x 1 diagonal block
            // Interchange rows K and IPIV(K).
            KP = IPIV( K )
            if (KP != K) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            K=K+1
         } else {
            // 2 x 2 diagonal block
            // Interchange rows K and -IPIV(K+1).
            KP = -IPIV( K+1 )
            IF( KP == -IPIV( K ) ) CALL ZSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )
            K=K+2
         }
        }

*  Compute (L \P**T * B) -> B    [ (L \P**T * B) ]

        ztrsm('L','L','N','U',N,NRHS,ONE,A,LDA,B,LDB);

*  Compute D \ B -> B   [ D \ (L \P**T * B) ]

         I=1
         DO WHILE ( I .LE. N )
            if ( IPIV(I) .GT. 0 ) {
              S = DBLE( ONE ) / DBLE( A( I, I ) )
              zdscal(NRHS, S, B( I, 1 ), LDB );
            } else {
                  AKM1K = WORK(I)
                  AKM1 = A( I, I ) / DCONJG( AKM1K )
                  AK = A( I+1, I+1 ) / AKM1K
                  DENOM = AKM1*AK - ONE
                  for (J = 1; J <= NRHS; J++) { // 25
                     BKM1 = B( I, J ) / DCONJG( AKM1K )
                     BK = B( I+1, J ) / AKM1K
                     B( I, J ) = ( AK*BKM1-BK ) / DENOM
                     B( I+1, J ) = ( AKM1*BK-BKM1 ) / DENOM
                 } // 25
                  I = I + 1
            }
            I = I + 1
         }

*  Compute (L**H \ B) -> B   [ L**H \ (D \ (L \P**T * B) ) ]

        ztrsm('L','L','C','U',N,NRHS,ONE,A,LDA,B,LDB);

        // P * B  [ P * (L**H \ (D \ (L \P**T * B) )) ]

        K=N
        DO WHILE ( K .GE. 1 )
         if ( IPIV( K ).GT.0 ) {
            // 1 x 1 diagonal block
            // Interchange rows K and IPIV(K).
            KP = IPIV( K )
            if (KP != K) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            K=K-1
         } else {
            // 2 x 2 diagonal block
            // Interchange rows K-1 and -IPIV(K).
            KP = -IPIV( K )
            IF( K.GT.1 && KP == -IPIV( K-1 ) ) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K=K-2
         }
        }

      }

      // Revert A

      zsyconv(UPLO, 'R', N, A, LDA, IPIV, WORK, IINFO );

      RETURN

      // End of ZHETRS2

      }
