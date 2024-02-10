      void slagsy(N, K, D, A, LDA, ISEED, WORK, Box<int> INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, K, LDA, N;
      int                ISEED( 4 );
      double               A( LDA, * ), D( * ), WORK( * );
      // ..

      double               ZERO, ONE, HALF;
      const              ZERO = 0.0, ONE = 1.0, HALF = 0.5 ;
      int                I, J;
      double               ALPHA, TAU, WA, WB, WN;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SGEMV, SGER, SLARNV, SSCAL, SSYMV, SSYR2, XERBLA
      // ..
      // .. External Functions ..
      //- REAL               SDOT, SNRM2;
      // EXTERNAL SDOT, SNRM2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SIGN

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
         xerbla('SLAGSY', -INFO );
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

      // Generate lower triangle of symmetric matrix

      for (I = N - 1; I >= 1; I--) { // 40

         // generate random reflection

         slarnv(3, ISEED, N-I+1, WORK );
         WN = SNRM2( N-I+1, WORK, 1 );
         WA = sign( WN, WORK( 1 ) );
         if ( WN == ZERO ) {
            TAU = ZERO;
         } else {
            WB = WORK( 1 ) + WA;
            sscal(N-I, ONE / WB, WORK( 2 ), 1 );
            WORK[1] = ONE;
            TAU = WB / WA;
         }

         // apply random reflection to A(i:n,i:n) from the left
         // and the right

         // compute  y := tau * A * u

         ssymv('Lower', N-I+1, TAU, A( I, I ), LDA, WORK, 1, ZERO, WORK( N+1 ), 1 );

         // compute  v := y - 1/2 * tau * ( y, u ) * u

         ALPHA = -HALF*TAU*SDOT( N-I+1, WORK( N+1 ), 1, WORK, 1 );
         saxpy(N-I+1, ALPHA, WORK, 1, WORK( N+1 ), 1 );

         // apply the transformation as a rank-2 update to A(i:n,i:n)

         ssyr2('Lower', N-I+1, -ONE, WORK, 1, WORK( N+1 ), 1, A( I, I ), LDA );
      } // 40

      // Reduce number of subdiagonals to K

      for (I = 1; I <= N - 1 - K; I++) { // 60

         // generate reflection to annihilate A(k+i+1:n,i)

         WN = SNRM2( N-K-I+1, A( K+I, I ), 1 );
         WA = sign( WN, A( K+I, I ) );
         if ( WN == ZERO ) {
            TAU = ZERO;
         } else {
            WB = A( K+I, I ) + WA;
            sscal(N-K-I, ONE / WB, A( K+I+1, I ), 1 );
            A[K+I][I] = ONE;
            TAU = WB / WA;
         }

         // apply reflection to A(k+i:n,i+1:k+i-1) from the left

         sgemv('Transpose', N-K-I+1, K-1, ONE, A( K+I, I+1 ), LDA, A( K+I, I ), 1, ZERO, WORK, 1 );
         sger(N-K-I+1, K-1, -TAU, A( K+I, I ), 1, WORK, 1, A( K+I, I+1 ), LDA );

         // apply reflection to A(k+i:n,k+i:n) from the left and the right

         // compute  y := tau * A * u

         ssymv('Lower', N-K-I+1, TAU, A( K+I, K+I ), LDA, A( K+I, I ), 1, ZERO, WORK, 1 );

         // compute  v := y - 1/2 * tau * ( y, u ) * u

         ALPHA = -HALF*TAU*SDOT( N-K-I+1, WORK, 1, A( K+I, I ), 1 );
         saxpy(N-K-I+1, ALPHA, A( K+I, I ), 1, WORK, 1 );

         // apply symmetric rank-2 update to A(k+i:n,k+i:n)

         ssyr2('Lower', N-K-I+1, -ONE, A( K+I, I ), 1, WORK, 1, A( K+I, K+I ), LDA );

         A[K+I][I] = -WA;
         for (J = K + I + 1; J <= N; J++) { // 50
            A[J][I] = ZERO;
         } // 50
      } // 60

      // Store full symmetric matrix

      for (J = 1; J <= N; J++) { // 80
         for (I = J + 1; I <= N; I++) { // 70
            A[J][I] = A( I, J );
         } // 70
      } // 80
      }
