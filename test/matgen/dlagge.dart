      void dlagge(M, N, KL, KU, D, final Matrix<double> A, final int LDA, final Array<int> ISEED, final Array<double> _WORK, final Box<int> INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, KL, KU, LDA, M, N;
      int                ISEED( 4 );
      double             A( LDA, * ), D( * ), WORK( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, J;
      double             TAU, WA, WB, WN;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMV, DGER, DLARNV, DSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SIGN
      // ..
      // .. External Functions ..
      //- double             DNRM2;
      // EXTERNAL DNRM2

      // Test the input arguments

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( KL < 0 || KL > M-1 ) {
         INFO = -3;
      } else if ( KU < 0 || KU > N-1 ) {
         INFO = -4;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -7;
      }
      if ( INFO < 0 ) {
         xerbla('DLAGGE', -INFO );
         return;
      }

      // initialize A to diagonal matrix

      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= M; I++) { // 10
            A[I][J] = ZERO;
         } // 10
      } // 20
      for (I = 1; I <= min( M, N ); I++) { // 30
         A[I][I] = D( I );
      } // 30

      // Quick exit if the user wants a diagonal matrix

      if(( KL == 0 ) && ( KU == 0)) return;

      // pre- and post-multiply A by random orthogonal matrices

      for (I = min( M, N ); I >= 1; I--) { // 40
         if ( I < M ) {

            // generate random reflection

            dlarnv(3, ISEED, M-I+1, WORK );
            WN = dnrm2( M-I+1, WORK, 1 );
            WA = sign( WN, WORK( 1 ) );
            if ( WN == ZERO ) {
               TAU = ZERO;
            } else {
               WB = WORK( 1 ) + WA;
               dscal(M-I, ONE / WB, WORK( 2 ), 1 );
               WORK[1] = ONE;
               TAU = WB / WA;
            }

            // multiply A(i:m,i:n) by random reflection from the left

            dgemv('Transpose', M-I+1, N-I+1, ONE, A( I, I ), LDA, WORK, 1, ZERO, WORK( M+1 ), 1 );
            dger(M-I+1, N-I+1, -TAU, WORK, 1, WORK( M+1 ), 1, A( I, I ), LDA );
         }
         if ( I < N ) {

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

            // multiply A(i:m,i:n) by random reflection from the right

            dgemv('No transpose', M-I+1, N-I+1, ONE, A( I, I ), LDA, WORK, 1, ZERO, WORK( N+1 ), 1 );
            dger(M-I+1, N-I+1, -TAU, WORK( N+1 ), 1, WORK, 1, A( I, I ), LDA );
         }
      } // 40

      // Reduce number of subdiagonals to KL and number of superdiagonals
      // to KU

      for (I = 1; I <= max( M-1-KL, N-1-KU ); I++) { // 70
         if ( KL <= KU ) {

            // annihilate subdiagonal elements first (necessary if KL = 0)

            if ( I <= min( M-1-KL, N ) ) {

               // generate reflection to annihilate A(kl+i+1:m,i)

               WN = dnrm2( M-KL-I+1, A( KL+I, I ), 1 );
               WA = sign( WN, A( KL+I, I ) );
               if ( WN == ZERO ) {
                  TAU = ZERO;
               } else {
                  WB = A( KL+I, I ) + WA;
                  dscal(M-KL-I, ONE / WB, A( KL+I+1, I ), 1 );
                  A[KL+I][I] = ONE;
                  TAU = WB / WA;
               }

               // apply reflection to A(kl+i:m,i+1:n) from the left

               dgemv('Transpose', M-KL-I+1, N-I, ONE, A( KL+I, I+1 ), LDA, A( KL+I, I ), 1, ZERO, WORK, 1 );
               dger(M-KL-I+1, N-I, -TAU, A( KL+I, I ), 1, WORK, 1, A( KL+I, I+1 ), LDA );
               A[KL+I][I] = -WA;
            }

            if ( I <= min( N-1-KU, M ) ) {

               // generate reflection to annihilate A(i,ku+i+1:n)

               WN = dnrm2( N-KU-I+1, A( I, KU+I ), LDA );
               WA = sign( WN, A( I, KU+I ) );
               if ( WN == ZERO ) {
                  TAU = ZERO;
               } else {
                  WB = A( I, KU+I ) + WA;
                  dscal(N-KU-I, ONE / WB, A( I, KU+I+1 ), LDA );
                  A[I][KU+I] = ONE;
                  TAU = WB / WA;
               }

               // apply reflection to A(i+1:m,ku+i:n) from the right

               dgemv('No transpose', M-I, N-KU-I+1, ONE, A( I+1, KU+I ), LDA, A( I, KU+I ), LDA, ZERO, WORK, 1 );
               dger(M-I, N-KU-I+1, -TAU, WORK, 1, A( I, KU+I ), LDA, A( I+1, KU+I ), LDA );
               A[I][KU+I] = -WA;
            }
         } else {

            // annihilate superdiagonal elements first (necessary if
            // KU = 0)

            if ( I <= min( N-1-KU, M ) ) {

               // generate reflection to annihilate A(i,ku+i+1:n)

               WN = dnrm2( N-KU-I+1, A( I, KU+I ), LDA );
               WA = sign( WN, A( I, KU+I ) );
               if ( WN == ZERO ) {
                  TAU = ZERO;
               } else {
                  WB = A( I, KU+I ) + WA;
                  dscal(N-KU-I, ONE / WB, A( I, KU+I+1 ), LDA );
                  A[I][KU+I] = ONE;
                  TAU = WB / WA;
               }

               // apply reflection to A(i+1:m,ku+i:n) from the right

               dgemv('No transpose', M-I, N-KU-I+1, ONE, A( I+1, KU+I ), LDA, A( I, KU+I ), LDA, ZERO, WORK, 1 );
               dger(M-I, N-KU-I+1, -TAU, WORK, 1, A( I, KU+I ), LDA, A( I+1, KU+I ), LDA );
               A[I][KU+I] = -WA;
            }

            if ( I <= min( M-1-KL, N ) ) {

               // generate reflection to annihilate A(kl+i+1:m,i)

               WN = dnrm2( M-KL-I+1, A( KL+I, I ), 1 );
               WA = sign( WN, A( KL+I, I ) );
               if ( WN == ZERO ) {
                  TAU = ZERO;
               } else {
                  WB = A( KL+I, I ) + WA;
                  dscal(M-KL-I, ONE / WB, A( KL+I+1, I ), 1 );
                  A[KL+I][I] = ONE;
                  TAU = WB / WA;
               }

               // apply reflection to A(kl+i:m,i+1:n) from the left

               dgemv('Transpose', M-KL-I+1, N-I, ONE, A( KL+I, I+1 ), LDA, A( KL+I, I ), 1, ZERO, WORK, 1 );
               dger(M-KL-I+1, N-I, -TAU, A( KL+I, I ), 1, WORK, 1, A( KL+I, I+1 ), LDA );
               A[KL+I][I] = -WA;
            }
         }

         if (I <= N) {
            for (J = KL + I + 1; J <= M; J++) { // 50
               A[J][I] = ZERO;
            } // 50
         }

         if (I <= M) {
            for (J = KU + I + 1; J <= N; J++) { // 60
               A[I][J] = ZERO;
            } // 60
         }
      } // 70
      }
