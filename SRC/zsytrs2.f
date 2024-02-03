      void zsytrs2(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      Complex       A( LDA, * ), B( LDB, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      Complex         ONE;
      const              ONE = (1.0,0.0) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, IINFO, J, K, KP;
      Complex         AK, AKM1, AKM1K, BK, BKM1, DENOM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZSCAL, ZSYCONV, ZSWAP, ZTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -8;
      }
      if ( INFO != 0 ) {
         xerbla('ZSYTRS2', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) return;

      // Convert A

      zsyconv(UPLO, 'C', N, A, LDA, IPIV, WORK, IINFO );

      if ( UPPER ) {

         // Solve A*X = B, where A = U*D*U**T.

        // P**T * B
        K=N;
        DO WHILE ( K >= 1 );
         if ( IPIV( K ) > 0 ) {
            // 1 x 1 diagonal block
            // Interchange rows K and IPIV(K).
            KP = IPIV( K );
            if (KP != K) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            K=K-1;
         } else {
            // 2 x 2 diagonal block
            // Interchange rows K-1 and -IPIV(K).
            KP = -IPIV( K );
            if( KP == -IPIV( K-1 ) ) CALL ZSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB );
            K=K-2;
         }
        }

// Compute (U \P**T * B) -> B    [ (U \P**T * B) ]

        ztrsm('L','U','N','U',N,NRHS,ONE,A,LDA,B,LDB);

// Compute D \ B -> B   [ D \ (U \P**T * B) ]

         I=N;
         DO WHILE ( I >= 1 );
            if ( IPIV(I) > 0 ) {
              zscal(NRHS, ONE / A( I, I ), B( I, 1 ), LDB );
            } else if ( I > 1) {
               if ( IPIV(I-1) == IPIV(I) ) {
                  AKM1K = WORK(I);
                  AKM1 = A( I-1, I-1 ) / AKM1K;
                  AK = A( I, I ) / AKM1K;
                  DENOM = AKM1*AK - ONE;
                  for (J = 1; J <= NRHS; J++) { // 15
                     BKM1 = B( I-1, J ) / AKM1K;
                     BK = B( I, J ) / AKM1K;
                     B( I-1, J ) = ( AK*BKM1-BK ) / DENOM;
                     B( I, J ) = ( AKM1*BK-BKM1 ) / DENOM;
                 } // 15
               I = I - 1;
               }
            }
            I = I - 1;
         }

       // Compute (U**T \ B) -> B   [ U**T \ (D \ (U \P**T * B) ) ]

         ztrsm('L','U','T','U',N,NRHS,ONE,A,LDA,B,LDB);

        // P * B  [ P * (U**T \ (D \ (U \P**T * B) )) ]

        K=1;
        DO WHILE ( K <= N );
         if ( IPIV( K ) > 0 ) {
            // 1 x 1 diagonal block
            // Interchange rows K and IPIV(K).
            KP = IPIV( K );
            if (KP != K) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            K=K+1;
         } else {
            // 2 x 2 diagonal block
            // Interchange rows K-1 and -IPIV(K).
            KP = -IPIV( K );
            if( K < N && KP == -IPIV( K+1 ) ) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            K=K+2;
         }
        }

      } else {

         // Solve A*X = B, where A = L*D*L**T.

        // P**T * B
        K=1;
        DO WHILE ( K <= N );
         if ( IPIV( K ) > 0 ) {
            // 1 x 1 diagonal block
            // Interchange rows K and IPIV(K).
            KP = IPIV( K );
            if (KP != K) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            K=K+1;
         } else {
            // 2 x 2 diagonal block
            // Interchange rows K and -IPIV(K+1).
            KP = -IPIV( K+1 );
            if( KP == -IPIV( K ) ) CALL ZSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB );
            K=K+2;
         }
        }

// Compute (L \P**T * B) -> B    [ (L \P**T * B) ]

        ztrsm('L','L','N','U',N,NRHS,ONE,A,LDA,B,LDB);

// Compute D \ B -> B   [ D \ (L \P**T * B) ]

         I=1;
         DO WHILE ( I <= N );
            if ( IPIV(I) > 0 ) {
              zscal(NRHS, ONE / A( I, I ), B( I, 1 ), LDB );
            } else {
                  AKM1K = WORK(I);
                  AKM1 = A( I, I ) / AKM1K;
                  AK = A( I+1, I+1 ) / AKM1K;
                  DENOM = AKM1*AK - ONE;
                  for (J = 1; J <= NRHS; J++) { // 25
                     BKM1 = B( I, J ) / AKM1K;
                     BK = B( I+1, J ) / AKM1K;
                     B( I, J ) = ( AK*BKM1-BK ) / DENOM;
                     B( I+1, J ) = ( AKM1*BK-BKM1 ) / DENOM;
                 } // 25
                  I = I + 1;
            }
            I = I + 1;
         }

// Compute (L**T \ B) -> B   [ L**T \ (D \ (L \P**T * B) ) ]

        ztrsm('L','L','T','U',N,NRHS,ONE,A,LDA,B,LDB);

        // P * B  [ P * (L**T \ (D \ (L \P**T * B) )) ]

        K=N;
        DO WHILE ( K >= 1 );
         if ( IPIV( K ) > 0 ) {
            // 1 x 1 diagonal block
            // Interchange rows K and IPIV(K).
            KP = IPIV( K );
            if (KP != K) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            K=K-1;
         } else {
            // 2 x 2 diagonal block
            // Interchange rows K-1 and -IPIV(K).
            KP = -IPIV( K );
            if( K > 1 && KP == -IPIV( K-1 ) ) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            K=K-2;
         }
        }

      }

      // Revert A

      zsyconv(UPLO, 'R', N, A, LDA, IPIV, WORK, IINFO );

      return;

      // End of ZSYTRS2

      }
