      void dlarge(N, A, LDA, ISEED, WORK, INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, N;
      int                ISEED( 4 );
      double             A( LDA, * ), WORK( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I;
      double             TAU, WA, WB, WN;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMV, DGER, DLARNV, DSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SIGN
      // ..
      // .. External Functions ..
      //- double             DNRM2;
      // EXTERNAL DNRM2

      // Test the input arguments

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -3;
      }
      if ( INFO < 0 ) {
         xerbla('DLARGE', -INFO );
         return;
      }

      // pre- and post-multiply A by random orthogonal matrix

      for (I = N; I >= 1; I--) { // 10

         // generate random reflection

         dlarnv(3, ISEED, N-I+1, WORK );
         WN = dnrm2( N-I+1, WORK, 1 );
         WA = sign( WN, WORK( 1 ) );
         if ( WN == ZERO ) {
            TAU = ZERO;
         } else {
            WB = WORK( 1 ) + WA;
            dscal(N-I, ONE / WB, WORK( 2 ), 1 );
            WORK[1] = ONE;
            TAU = WB / WA;
         }

         // multiply A(i:n,1:n) by random reflection from the left

         dgemv('Transpose', N-I+1, N, ONE, A( I, 1 ), LDA, WORK, 1, ZERO, WORK( N+1 ), 1 );
         dger(N-I+1, N, -TAU, WORK, 1, WORK( N+1 ), 1, A( I, 1 ), LDA );

         // multiply A(1:n,i:n) by random reflection from the right

         dgemv('No transpose', N, N-I+1, ONE, A( 1, I ), LDA, WORK, 1, ZERO, WORK( N+1 ), 1 );
         dger(N, N-I+1, -TAU, WORK( N+1 ), 1, WORK, 1, A( 1, I ), LDA );
      } // 10
      }
