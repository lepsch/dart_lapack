      void chegs2(ITYPE, UPLO, N, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, ITYPE, LDA, LDB, N;
      Complex            A( LDA, * ), B( LDB, * );
      // ..

      double               ONE, HALF;
      const              ONE = 1.0, HALF = 0.5 ;
      Complex            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      bool               UPPER;
      int                K;
      double               AKK, BKK;
      Complex            CT;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CHER2, CLACGV, CSSCAL, CTRMV, CTRSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( ITYPE < 1 || ITYPE > 3 ) {
         INFO = -1;
      } else if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('CHEGS2', -INFO );
         return;
      }

      if ( ITYPE == 1 ) {
         if ( UPPER ) {

            // Compute inv(U**H)*A*inv(U)

            for (K = 1; K <= N; K++) { // 10

               // Update the upper triangle of A(k:n,k:n)

               AKK = double( A( K, K ) );
               BKK = double( B( K, K ) );
               AKK = AKK / BKK**2;
               A[K][K] = AKK;
               if ( K < N ) {
                  csscal(N-K, ONE / BKK, A( K, K+1 ), LDA );
                  CT = -HALF*AKK;
                  clacgv(N-K, A( K, K+1 ), LDA );
                  clacgv(N-K, B( K, K+1 ), LDB );
                  caxpy(N-K, CT, B( K, K+1 ), LDB, A( K, K+1 ), LDA );
                  cher2(UPLO, N-K, -CONE, A( K, K+1 ), LDA, B( K, K+1 ), LDB, A( K+1, K+1 ), LDA );
                  caxpy(N-K, CT, B( K, K+1 ), LDB, A( K, K+1 ), LDA );
                  clacgv(N-K, B( K, K+1 ), LDB );
                  ctrsv(UPLO, 'Conjugate transpose', 'Non-unit', N-K, B( K+1, K+1 ), LDB, A( K, K+1 ), LDA );
                  clacgv(N-K, A( K, K+1 ), LDA );
               }
            } // 10
         } else {

            // Compute inv(L)*A*inv(L**H)

            for (K = 1; K <= N; K++) { // 20

               // Update the lower triangle of A(k:n,k:n)

               AKK = double( A( K, K ) );
               BKK = double( B( K, K ) );
               AKK = AKK / BKK**2;
               A[K][K] = AKK;
               if ( K < N ) {
                  csscal(N-K, ONE / BKK, A( K+1, K ), 1 );
                  CT = -HALF*AKK;
                  caxpy(N-K, CT, B( K+1, K ), 1, A( K+1, K ), 1 );
                  cher2(UPLO, N-K, -CONE, A( K+1, K ), 1, B( K+1, K ), 1, A( K+1, K+1 ), LDA );
                  caxpy(N-K, CT, B( K+1, K ), 1, A( K+1, K ), 1 );
                  ctrsv(UPLO, 'No transpose', 'Non-unit', N-K, B( K+1, K+1 ), LDB, A( K+1, K ), 1 );
               }
            } // 20
         }
      } else {
         if ( UPPER ) {

            // Compute U*A*U**H

            for (K = 1; K <= N; K++) { // 30

               // Update the upper triangle of A(1:k,1:k)

               AKK = double( A( K, K ) );
               BKK = double( B( K, K ) );
               ctrmv(UPLO, 'No transpose', 'Non-unit', K-1, B, LDB, A( 1, K ), 1 );
               CT = HALF*AKK;
               caxpy(K-1, CT, B( 1, K ), 1, A( 1, K ), 1 );
               cher2(UPLO, K-1, CONE, A( 1, K ), 1, B( 1, K ), 1, A, LDA );
               caxpy(K-1, CT, B( 1, K ), 1, A( 1, K ), 1 );
               csscal(K-1, BKK, A( 1, K ), 1 );
               A[K][K] = AKK*BKK**2;
            } // 30
         } else {

            // Compute L**H *A*L

            for (K = 1; K <= N; K++) { // 40

               // Update the lower triangle of A(1:k,1:k)

               AKK = double( A( K, K ) );
               BKK = double( B( K, K ) );
               clacgv(K-1, A( K, 1 ), LDA );
               ctrmv(UPLO, 'Conjugate transpose', 'Non-unit', K-1, B, LDB, A( K, 1 ), LDA );
               CT = HALF*AKK;
               clacgv(K-1, B( K, 1 ), LDB );
               caxpy(K-1, CT, B( K, 1 ), LDB, A( K, 1 ), LDA );
               cher2(UPLO, K-1, CONE, A( K, 1 ), LDA, B( K, 1 ), LDB, A, LDA );
               caxpy(K-1, CT, B( K, 1 ), LDB, A( K, 1 ), LDA );
               clacgv(K-1, B( K, 1 ), LDB );
               csscal(K-1, BKK, A( K, 1 ), LDA );
               clacgv(K-1, A( K, 1 ), LDA );
               A[K][K] = AKK*BKK**2;
            } // 40
         }
      }
      }
