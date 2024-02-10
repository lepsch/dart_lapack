      void ssytrs_rook(UPLO, N, NRHS, final Matrix<double> A, final int LDA, IPIV, final Matrix<double> B, final int LDB, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, LDB, N, NRHS;
      int                IPIV( * );
      double               A( LDA, * ), B( LDB, * );
      // ..

      double               ONE;
      const              ONE = 1.0 ;
      bool               UPPER;
      int                J, K, KP;
      double               AK, AKM1, AKM1K, BK, BKM1, DENOM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMV, SGER, SSCAL, SSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
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
         xerbla('SSYTRS_ROOK', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) return;

      if ( UPPER ) {

         // Solve A*X = B, where A = U*D*U**T.

         // First solve U*D*X = B, overwriting B with X.

         // K is the main loop index, decreasing from N to 1 in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = N;
         } // 10

         // If K < 1, exit from loop.

         if (K < 1) GO TO 30;

         if ( IPIV( K ) > 0 ) {

            // 1 x 1 diagonal block

            // Interchange rows K and IPIV(K).

            KP = IPIV( K );
            if (KP != K) sswap( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );

            // Multiply by inv(U(K)), where U(K) is the transformation
            // stored in column K of A.

            sger(K-1, NRHS, -ONE, A( 1, K ), 1, B( K, 1 ), LDB, B( 1, 1 ), LDB );

            // Multiply by the inverse of the diagonal block.

            sscal(NRHS, ONE / A( K, K ), B( K, 1 ), LDB );
            K = K - 1;
         } else {

            // 2 x 2 diagonal block

            // Interchange rows K and -IPIV(K) THEN K-1 and -IPIV(K-1)

            KP = -IPIV( K );
            if (KP != K) sswap( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );

            KP = -IPIV( K-1 );
            if (KP != K-1) sswap( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB );

            // Multiply by inv(U(K)), where U(K) is the transformation
            // stored in columns K-1 and K of A.

            if ( K > 2 ) {
               sger(K-2, NRHS, -ONE, A( 1, K ), 1, B( K, 1 ), LDB, B( 1, 1 ), LDB );
               sger(K-2, NRHS, -ONE, A( 1, K-1 ), 1, B( K-1, 1 ), LDB, B( 1, 1 ), LDB );
            }

            // Multiply by the inverse of the diagonal block.

            AKM1K = A( K-1, K );
            AKM1 = A( K-1, K-1 ) / AKM1K;
            AK = A( K, K ) / AKM1K;
            DENOM = AKM1*AK - ONE;
            for (J = 1; J <= NRHS; J++) { // 20
               BKM1 = B( K-1, J ) / AKM1K;
               BK = B( K, J ) / AKM1K;
               B[K-1][J] = ( AK*BKM1-BK ) / DENOM;
               B[K][J] = ( AKM1*BK-BKM1 ) / DENOM;
            } // 20
            K = K - 2;
         }

         GO TO 10;
         } // 30

         // Next solve U**T *X = B, overwriting B with X.

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = 1;
         } // 40

         // If K > N, exit from loop.

         if (K > N) GO TO 50;

         if ( IPIV( K ) > 0 ) {

            // 1 x 1 diagonal block

            // Multiply by inv(U**T(K)), where U(K) is the transformation
            // stored in column K of A.

            if (K > 1) sgemv( 'Transpose', K-1, NRHS, -ONE, B, LDB, A( 1, K ), 1, ONE, B( K, 1 ), LDB );

            // Interchange rows K and IPIV(K).

            KP = IPIV( K );
            if (KP != K) sswap( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            K = K + 1;
         } else {

            // 2 x 2 diagonal block

            // Multiply by inv(U**T(K+1)), where U(K+1) is the transformation
            // stored in columns K and K+1 of A.

            if ( K > 1 ) {
               sgemv('Transpose', K-1, NRHS, -ONE, B, LDB, A( 1, K ), 1, ONE, B( K, 1 ), LDB );
               sgemv('Transpose', K-1, NRHS, -ONE, B, LDB, A( 1, K+1 ), 1, ONE, B( K+1, 1 ), LDB );
            }

            // Interchange rows K and -IPIV(K) THEN K+1 and -IPIV(K+1).

            KP = -IPIV( K );
            if (KP != K) sswap( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );

            KP = -IPIV( K+1 );
            if (KP != K+1) sswap( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB );

            K = K + 2;
         }

         GO TO 40;
         } // 50

      } else {

         // Solve A*X = B, where A = L*D*L**T.

         // First solve L*D*X = B, overwriting B with X.

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = 1;
         } // 60

         // If K > N, exit from loop.

         if (K > N) GO TO 80;

         if ( IPIV( K ) > 0 ) {

            // 1 x 1 diagonal block

            // Interchange rows K and IPIV(K).

            KP = IPIV( K );
            if (KP != K) sswap( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );

            // Multiply by inv(L(K)), where L(K) is the transformation
            // stored in column K of A.

            if (K < N) sger( N-K, NRHS, -ONE, A( K+1, K ), 1, B( K, 1 ), LDB, B( K+1, 1 ), LDB );

            // Multiply by the inverse of the diagonal block.

            sscal(NRHS, ONE / A( K, K ), B( K, 1 ), LDB );
            K = K + 1;
         } else {

            // 2 x 2 diagonal block

            // Interchange rows K and -IPIV(K) THEN K+1 and -IPIV(K+1)

            KP = -IPIV( K );
            if (KP != K) sswap( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );

            KP = -IPIV( K+1 );
            if (KP != K+1) sswap( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB );

            // Multiply by inv(L(K)), where L(K) is the transformation
            // stored in columns K and K+1 of A.

            if ( K < N-1 ) {
               sger(N-K-1, NRHS, -ONE, A( K+2, K ), 1, B( K, 1 ), LDB, B( K+2, 1 ), LDB );
               sger(N-K-1, NRHS, -ONE, A( K+2, K+1 ), 1, B( K+1, 1 ), LDB, B( K+2, 1 ), LDB );
            }

            // Multiply by the inverse of the diagonal block.

            AKM1K = A( K+1, K );
            AKM1 = A( K, K ) / AKM1K;
            AK = A( K+1, K+1 ) / AKM1K;
            DENOM = AKM1*AK - ONE;
            for (J = 1; J <= NRHS; J++) { // 70
               BKM1 = B( K, J ) / AKM1K;
               BK = B( K+1, J ) / AKM1K;
               B[K][J] = ( AK*BKM1-BK ) / DENOM;
               B[K+1][J] = ( AKM1*BK-BKM1 ) / DENOM;
            } // 70
            K = K + 2;
         }

         GO TO 60;
         } // 80

         // Next solve L**T *X = B, overwriting B with X.

         // K is the main loop index, decreasing from N to 1 in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = N;
         } // 90

         // If K < 1, exit from loop.

         if (K < 1) GO TO 100;

         if ( IPIV( K ) > 0 ) {

            // 1 x 1 diagonal block

            // Multiply by inv(L**T(K)), where L(K) is the transformation
            // stored in column K of A.

            if (K < N) sgemv( 'Transpose', N-K, NRHS, -ONE, B( K+1, 1 ), LDB, A( K+1, K ), 1, ONE, B( K, 1 ), LDB );

            // Interchange rows K and IPIV(K).

            KP = IPIV( K );
            if (KP != K) sswap( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            K = K - 1;
         } else {

            // 2 x 2 diagonal block

            // Multiply by inv(L**T(K-1)), where L(K-1) is the transformation
            // stored in columns K-1 and K of A.

            if ( K < N ) {
               sgemv('Transpose', N-K, NRHS, -ONE, B( K+1, 1 ), LDB, A( K+1, K ), 1, ONE, B( K, 1 ), LDB );
               sgemv('Transpose', N-K, NRHS, -ONE, B( K+1, 1 ), LDB, A( K+1, K-1 ), 1, ONE, B( K-1, 1 ), LDB );
            }

            // Interchange rows K and -IPIV(K) THEN K-1 and -IPIV(K-1)

            KP = -IPIV( K );
            if (KP != K) sswap( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );

            KP = -IPIV( K-1 );
            if (KP != K-1) sswap( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB );

            K = K - 2;
         }

         GO TO 90;
         } // 100
      }

      }
