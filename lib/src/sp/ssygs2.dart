      void ssygs2(final int ITYPE, final int UPLO, final int N, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, ITYPE, LDA, LDB, N;
      double               A( LDA, * ), B( LDB, * );
      // ..

      double               ONE, HALF;
      const              ONE = 1.0, HALF = 0.5 ;
      bool               UPPER;
      int                K;
      double               AKK, BKK, CT;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SSCAL, SSYR2, STRMV, STRSV, XERBLA
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
         xerbla('SSYGS2', -INFO );
         return;
      }

      if ( ITYPE == 1 ) {
         if ( UPPER ) {

            // Compute inv(U**T)*A*inv(U)

            for (K = 1; K <= N; K++) { // 10

               // Update the upper triangle of A(k:n,k:n)

               AKK = A( K, K );
               BKK = B( K, K );
               AKK = AKK / BKK**2;
               A[K][K] = AKK;
               if ( K < N ) {
                  sscal(N-K, ONE / BKK, A( K, K+1 ), LDA );
                  CT = -HALF*AKK;
                  saxpy(N-K, CT, B( K, K+1 ), LDB, A( K, K+1 ), LDA );
                  ssyr2(UPLO, N-K, -ONE, A( K, K+1 ), LDA, B( K, K+1 ), LDB, A( K+1, K+1 ), LDA );
                  saxpy(N-K, CT, B( K, K+1 ), LDB, A( K, K+1 ), LDA );
                  strsv(UPLO, 'Transpose', 'Non-unit', N-K, B( K+1, K+1 ), LDB, A( K, K+1 ), LDA );
               }
            } // 10
         } else {

            // Compute inv(L)*A*inv(L**T)

            for (K = 1; K <= N; K++) { // 20

               // Update the lower triangle of A(k:n,k:n)

               AKK = A( K, K );
               BKK = B( K, K );
               AKK = AKK / BKK**2;
               A[K][K] = AKK;
               if ( K < N ) {
                  sscal(N-K, ONE / BKK, A( K+1, K ), 1 );
                  CT = -HALF*AKK;
                  saxpy(N-K, CT, B( K+1, K ), 1, A( K+1, K ), 1 );
                  ssyr2(UPLO, N-K, -ONE, A( K+1, K ), 1, B( K+1, K ), 1, A( K+1, K+1 ), LDA );
                  saxpy(N-K, CT, B( K+1, K ), 1, A( K+1, K ), 1 );
                  strsv(UPLO, 'No transpose', 'Non-unit', N-K, B( K+1, K+1 ), LDB, A( K+1, K ), 1 );
               }
            } // 20
         }
      } else {
         if ( UPPER ) {

            // Compute U*A*U**T

            for (K = 1; K <= N; K++) { // 30

               // Update the upper triangle of A(1:k,1:k)

               AKK = A( K, K );
               BKK = B( K, K );
               strmv(UPLO, 'No transpose', 'Non-unit', K-1, B, LDB, A( 1, K ), 1 );
               CT = HALF*AKK;
               saxpy(K-1, CT, B( 1, K ), 1, A( 1, K ), 1 );
               ssyr2(UPLO, K-1, ONE, A( 1, K ), 1, B( 1, K ), 1, A, LDA );
               saxpy(K-1, CT, B( 1, K ), 1, A( 1, K ), 1 );
               sscal(K-1, BKK, A( 1, K ), 1 );
               A[K][K] = AKK*BKK**2;
            } // 30
         } else {

            // Compute L**T *A*L

            for (K = 1; K <= N; K++) { // 40

               // Update the lower triangle of A(1:k,1:k)

               AKK = A( K, K );
               BKK = B( K, K );
               strmv(UPLO, 'Transpose', 'Non-unit', K-1, B, LDB, A( K, 1 ), LDA );
               CT = HALF*AKK;
               saxpy(K-1, CT, B( K, 1 ), LDB, A( K, 1 ), LDA );
               ssyr2(UPLO, K-1, ONE, A( K, 1 ), LDA, B( K, 1 ), LDB, A, LDA );
               saxpy(K-1, CT, B( K, 1 ), LDB, A( K, 1 ), LDA );
               sscal(K-1, BKK, A( K, 1 ), LDA );
               A[K][K] = AKK*BKK**2;
            } // 40
         }
      }
      }
