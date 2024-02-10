      void zlagsy(N, K, D, A, LDA, ISEED, WORK, INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, K, LDA, N;
      int                ISEED( 4 );
      double             D( * );
      Complex         A( LDA, * ), WORK( * );
      // ..

      Complex         ZERO, ONE, HALF;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ), HALF = ( 0.5, 0.0 ) ;
      int                I, II, J, JJ;
      double             WN;
      Complex         ALPHA, TAU, WA, WB;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZAXPY, ZGEMV, ZGERC, ZLACGV, ZLARNV, ZSCAL, ZSYMV
      // ..
      // .. External Functions ..
      //- double             DZNRM2;
      //- Complex         ZDOTC;
      // EXTERNAL DZNRM2, ZDOTC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX

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
         xerbla('ZLAGSY', -INFO );
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

      for (I = N - 1; I >= 1; I--) { // 60

         // generate random reflection

         zlarnv(3, ISEED, N-I+1, WORK );
         WN = DZNRM2( N-I+1, WORK, 1 );
         WA = ( WN / ( WORK( 1 ) ).abs() )*WORK( 1 );
         if ( WN == ZERO ) {
            TAU = ZERO;
         } else {
            WB = WORK( 1 ) + WA;
            zscal(N-I, ONE / WB, WORK( 2 ), 1 );
            WORK[1] = ONE;
            TAU = (WB / WA).toDouble();
         }

         // apply random reflection to A(i:n,i:n) from the left
         // and the right

         // compute  y := tau * A * conjg(u)

         zlacgv(N-I+1, WORK, 1 );
         zsymv('Lower', N-I+1, TAU, A( I, I ), LDA, WORK, 1, ZERO, WORK( N+1 ), 1 );
         zlacgv(N-I+1, WORK, 1 );

         // compute  v := y - 1/2 * tau * ( u, y ) * u

         ALPHA = -HALF*TAU*ZDOTC( N-I+1, WORK, 1, WORK( N+1 ), 1 );
         zaxpy(N-I+1, ALPHA, WORK, 1, WORK( N+1 ), 1 );

         // apply the transformation as a rank-2 update to A(i:n,i:n)

         // CALL ZSYR2( 'Lower', N-I+1, -ONE, WORK, 1, WORK( N+1 ), 1,
         // $               A( I, I ), LDA )

         for (JJ = I; JJ <= N; JJ++) { // 50
            for (II = JJ; II <= N; II++) { // 40
               A[II][JJ] = A( II, JJ ) - WORK( II-I+1 )*WORK( N+JJ-I+1 ) - WORK( N+II-I+1 )*WORK( JJ-I+1 );
            } // 40
         } // 50
      } // 60

      // Reduce number of subdiagonals to K

      for (I = 1; I <= N - 1 - K; I++) { // 100

         // generate reflection to annihilate A(k+i+1:n,i)

         WN = DZNRM2( N-K-I+1, A( K+I, I ), 1 );
         WA = ( WN / ( A( K+I, I ) ).abs() )*A( K+I, I );
         if ( WN == ZERO ) {
            TAU = ZERO;
         } else {
            WB = A( K+I, I ) + WA;
            zscal(N-K-I, ONE / WB, A( K+I+1, I ), 1 );
            A[K+I][I] = ONE;
            TAU = (WB / WA).toDouble();
         }

         // apply reflection to A(k+i:n,i+1:k+i-1) from the left

         zgemv('Conjugate transpose', N-K-I+1, K-1, ONE, A( K+I, I+1 ), LDA, A( K+I, I ), 1, ZERO, WORK, 1 );
         zgerc(N-K-I+1, K-1, -TAU, A( K+I, I ), 1, WORK, 1, A( K+I, I+1 ), LDA );

         // apply reflection to A(k+i:n,k+i:n) from the left and the right

         // compute  y := tau * A * conjg(u)

         zlacgv(N-K-I+1, A( K+I, I ), 1 );
         zsymv('Lower', N-K-I+1, TAU, A( K+I, K+I ), LDA, A( K+I, I ), 1, ZERO, WORK, 1 );
         zlacgv(N-K-I+1, A( K+I, I ), 1 );

         // compute  v := y - 1/2 * tau * ( u, y ) * u

         ALPHA = -HALF*TAU*ZDOTC( N-K-I+1, A( K+I, I ), 1, WORK, 1 );
         zaxpy(N-K-I+1, ALPHA, A( K+I, I ), 1, WORK, 1 );

         // apply symmetric rank-2 update to A(k+i:n,k+i:n)

         // CALL ZSYR2( 'Lower', N-K-I+1, -ONE, A( K+I, I ), 1, WORK, 1,
         // $               A( K+I, K+I ), LDA )

         for (JJ = K + I; JJ <= N; JJ++) { // 80
            for (II = JJ; II <= N; II++) { // 70
               A[II][JJ] = A( II, JJ ) - A( II, I )*WORK( JJ-K-I+1 ) - WORK( II-K-I+1 )*A( JJ, I );
            } // 70
         } // 80

         A[K+I][I] = -WA;
         for (J = K + I + 1; J <= N; J++) { // 90
            A[J][I] = ZERO;
         } // 90
      } // 100

      // Store full symmetric matrix

      for (J = 1; J <= N; J++) { // 120
         for (I = J + 1; I <= N; I++) { // 110
            A[J][I] = A( I, J );
         } // 110
      } // 120
      }
