      void claghe(N, K, D, final Matrix<double> A, final int LDA, final Array<int> ISEED, final Array<double> _WORK, final Box<int> INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, K, LDA, N;
      int                ISEED( 4 );
      double               D( * );
      Complex            A( LDA, * ), WORK( * );
      // ..

      Complex            ZERO, ONE, HALF;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ), HALF = ( 0.5, 0.0 ) ;
      int                I, J;
      double               WN;
      Complex            ALPHA, TAU, WA, WB;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CGEMV, CGERC, CHEMV, CHER2, CLARNV, CSCAL, XERBLA
      // ..
      // .. External Functions ..
      //- REAL               SCNRM2;
      //- COMPLEX            CDOTC;
      // EXTERNAL SCNRM2, CDOTC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, MAX, REAL

      // Test the input arguments

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
      } else if ( K < 0 || K > N-1 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      }
      if ( INFO < 0 ) {
         xerbla('CLAGHE', -INFO );
         return;
      }

      // initialize lower triangle of A to diagonal matrix

      for (J = 1; J <= N; J++) { // 20
         for (I = J + 1; I <= N; I++) { // 10
            A[I][J] = ZERO;
         } // 10
      } // 20
      for (I = 1; I <= N; I++) { // 30
         A[I][I] = D( I );
      } // 30

      // Generate lower triangle of hermitian matrix

      for (I = N - 1; I >= 1; I--) { // 40

         // generate random reflection

         clarnv(3, ISEED, N-I+1, WORK );
         WN = SCNRM2( N-I+1, WORK, 1 );
         WA = ( WN / ( WORK( 1 ) ).abs() )*WORK( 1 );
         if ( WN == ZERO ) {
            TAU = ZERO;
         } else {
            WB = WORK( 1 ) + WA;
            cscal(N-I, ONE / WB, WORK( 2 ), 1 );
            WORK[1] = ONE;
            TAU = double( WB / WA );
         }

         // apply random reflection to A(i:n,i:n) from the left
         // and the right

         // compute  y := tau * A * u

         chemv('Lower', N-I+1, TAU, A( I, I ), LDA, WORK, 1, ZERO, WORK( N+1 ), 1 );

         // compute  v := y - 1/2 * tau * ( y, u ) * u

         ALPHA = -HALF*TAU*CDOTC( N-I+1, WORK( N+1 ), 1, WORK, 1 );
         caxpy(N-I+1, ALPHA, WORK, 1, WORK( N+1 ), 1 );

         // apply the transformation as a rank-2 update to A(i:n,i:n)

         cher2('Lower', N-I+1, -ONE, WORK, 1, WORK( N+1 ), 1, A( I, I ), LDA );
      } // 40

      // Reduce number of subdiagonals to K

      for (I = 1; I <= N - 1 - K; I++) { // 60

         // generate reflection to annihilate A(k+i+1:n,i)

         WN = SCNRM2( N-K-I+1, A( K+I, I ), 1 );
         WA = ( WN / ( A( K+I, I ) ).abs() )*A( K+I, I );
         if ( WN == ZERO ) {
            TAU = ZERO;
         } else {
            WB = A( K+I, I ) + WA;
            cscal(N-K-I, ONE / WB, A( K+I+1, I ), 1 );
            A[K+I][I] = ONE;
            TAU = double( WB / WA );
         }

         // apply reflection to A(k+i:n,i+1:k+i-1) from the left

         cgemv('Conjugate transpose', N-K-I+1, K-1, ONE, A( K+I, I+1 ), LDA, A( K+I, I ), 1, ZERO, WORK, 1 );
         cgerc(N-K-I+1, K-1, -TAU, A( K+I, I ), 1, WORK, 1, A( K+I, I+1 ), LDA );

         // apply reflection to A(k+i:n,k+i:n) from the left and the right

         // compute  y := tau * A * u

         chemv('Lower', N-K-I+1, TAU, A( K+I, K+I ), LDA, A( K+I, I ), 1, ZERO, WORK, 1 );

         // compute  v := y - 1/2 * tau * ( y, u ) * u

         ALPHA = -HALF*TAU*CDOTC( N-K-I+1, WORK, 1, A( K+I, I ), 1 );
         caxpy(N-K-I+1, ALPHA, A( K+I, I ), 1, WORK, 1 );

         // apply hermitian rank-2 update to A(k+i:n,k+i:n)

         cher2('Lower', N-K-I+1, -ONE, A( K+I, I ), 1, WORK, 1, A( K+I, K+I ), LDA );

         A[K+I][I] = -WA;
         for (J = K + I + 1; J <= N; J++) { // 50
            A[J][I] = ZERO;
         } // 50
      } // 60

      // Store full hermitian matrix

      for (J = 1; J <= N; J++) { // 80
         for (I = J + 1; I <= N; I++) { // 70
            A[J][I] = CONJG( A( I, J ) );
         } // 70
      } // 80
      }
