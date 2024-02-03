      void slagge(M, N, KL, KU, D, A, LDA, ISEED, WORK, INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, KL, KU, LDA, M, N;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      REAL               A( LDA, * ), D( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               TAU, WA, WB, WN;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMV, SGER, SLARNV, SSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SIGN
      // ..
      // .. External Functions ..
      REAL               SNRM2;
      // EXTERNAL SNRM2
      // ..
      // .. Executable Statements ..

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
         xerbla('SLAGGE', -INFO );
         return;
      }

      // initialize A to diagonal matrix

      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= M; I++) { // 10
            A( I, J ) = ZERO;
         } // 10
      } // 20
      DO 30 I = 1, min( M, N );
         A( I, I ) = D( I );
      } // 30

      // Quick exit if the user wants a diagonal matrix

      if(( KL == 0 ) && ( KU == 0)) return;

      // pre- and post-multiply A by random orthogonal matrices

      DO 40 I = min( M, N ), 1, -1;
         if ( I < M ) {

            // generate random reflection

            slarnv(3, ISEED, M-I+1, WORK );
            WN = SNRM2( M-I+1, WORK, 1 );
            WA = SIGN( WN, WORK( 1 ) );
            if ( WN == ZERO ) {
               TAU = ZERO;
            } else {
               WB = WORK( 1 ) + WA;
               sscal(M-I, ONE / WB, WORK( 2 ), 1 );
               WORK( 1 ) = ONE;
               TAU = WB / WA;
            }

            // multiply A(i:m,i:n) by random reflection from the left

            sgemv('Transpose', M-I+1, N-I+1, ONE, A( I, I ), LDA, WORK, 1, ZERO, WORK( M+1 ), 1 );
            sger(M-I+1, N-I+1, -TAU, WORK, 1, WORK( M+1 ), 1, A( I, I ), LDA );
         }
         if ( I < N ) {

            // generate random reflection

            slarnv(3, ISEED, N-I+1, WORK );
            WN = SNRM2( N-I+1, WORK, 1 );
            WA = SIGN( WN, WORK( 1 ) );
            if ( WN == ZERO ) {
               TAU = ZERO;
            } else {
               WB = WORK( 1 ) + WA;
               sscal(N-I, ONE / WB, WORK( 2 ), 1 );
               WORK( 1 ) = ONE;
               TAU = WB / WA;
            }

            // multiply A(i:m,i:n) by random reflection from the right

            sgemv('No transpose', M-I+1, N-I+1, ONE, A( I, I ), LDA, WORK, 1, ZERO, WORK( N+1 ), 1 );
            sger(M-I+1, N-I+1, -TAU, WORK( N+1 ), 1, WORK, 1, A( I, I ), LDA );
         }
      } // 40

      // Reduce number of subdiagonals to KL and number of superdiagonals
      // to KU

      DO 70 I = 1, max( M-1-KL, N-1-KU );
         if ( KL <= KU ) {

            // annihilate subdiagonal elements first (necessary if KL = 0)

            if ( I <= min( M-1-KL, N ) ) {

               // generate reflection to annihilate A(kl+i+1:m,i)

               WN = SNRM2( M-KL-I+1, A( KL+I, I ), 1 );
               WA = SIGN( WN, A( KL+I, I ) );
               if ( WN == ZERO ) {
                  TAU = ZERO;
               } else {
                  WB = A( KL+I, I ) + WA;
                  sscal(M-KL-I, ONE / WB, A( KL+I+1, I ), 1 );
                  A( KL+I, I ) = ONE;
                  TAU = WB / WA;
               }

               // apply reflection to A(kl+i:m,i+1:n) from the left

               sgemv('Transpose', M-KL-I+1, N-I, ONE, A( KL+I, I+1 ), LDA, A( KL+I, I ), 1, ZERO, WORK, 1 );
               sger(M-KL-I+1, N-I, -TAU, A( KL+I, I ), 1, WORK, 1, A( KL+I, I+1 ), LDA );
               A( KL+I, I ) = -WA;
            }

            if ( I <= min( N-1-KU, M ) ) {

               // generate reflection to annihilate A(i,ku+i+1:n)

               WN = SNRM2( N-KU-I+1, A( I, KU+I ), LDA );
               WA = SIGN( WN, A( I, KU+I ) );
               if ( WN == ZERO ) {
                  TAU = ZERO;
               } else {
                  WB = A( I, KU+I ) + WA;
                  sscal(N-KU-I, ONE / WB, A( I, KU+I+1 ), LDA );
                  A( I, KU+I ) = ONE;
                  TAU = WB / WA;
               }

               // apply reflection to A(i+1:m,ku+i:n) from the right

               sgemv('No transpose', M-I, N-KU-I+1, ONE, A( I+1, KU+I ), LDA, A( I, KU+I ), LDA, ZERO, WORK, 1 );
               sger(M-I, N-KU-I+1, -TAU, WORK, 1, A( I, KU+I ), LDA, A( I+1, KU+I ), LDA );
               A( I, KU+I ) = -WA;
            }
         } else {

            // annihilate superdiagonal elements first (necessary if
            // KU = 0)

            if ( I <= min( N-1-KU, M ) ) {

               // generate reflection to annihilate A(i,ku+i+1:n)

               WN = SNRM2( N-KU-I+1, A( I, KU+I ), LDA );
               WA = SIGN( WN, A( I, KU+I ) );
               if ( WN == ZERO ) {
                  TAU = ZERO;
               } else {
                  WB = A( I, KU+I ) + WA;
                  sscal(N-KU-I, ONE / WB, A( I, KU+I+1 ), LDA );
                  A( I, KU+I ) = ONE;
                  TAU = WB / WA;
               }

               // apply reflection to A(i+1:m,ku+i:n) from the right

               sgemv('No transpose', M-I, N-KU-I+1, ONE, A( I+1, KU+I ), LDA, A( I, KU+I ), LDA, ZERO, WORK, 1 );
               sger(M-I, N-KU-I+1, -TAU, WORK, 1, A( I, KU+I ), LDA, A( I+1, KU+I ), LDA );
               A( I, KU+I ) = -WA;
            }

            if ( I <= min( M-1-KL, N ) ) {

               // generate reflection to annihilate A(kl+i+1:m,i)

               WN = SNRM2( M-KL-I+1, A( KL+I, I ), 1 );
               WA = SIGN( WN, A( KL+I, I ) );
               if ( WN == ZERO ) {
                  TAU = ZERO;
               } else {
                  WB = A( KL+I, I ) + WA;
                  sscal(M-KL-I, ONE / WB, A( KL+I+1, I ), 1 );
                  A( KL+I, I ) = ONE;
                  TAU = WB / WA;
               }

               // apply reflection to A(kl+i:m,i+1:n) from the left

               sgemv('Transpose', M-KL-I+1, N-I, ONE, A( KL+I, I+1 ), LDA, A( KL+I, I ), 1, ZERO, WORK, 1 );
               sger(M-KL-I+1, N-I, -TAU, A( KL+I, I ), 1, WORK, 1, A( KL+I, I+1 ), LDA );
               A( KL+I, I ) = -WA;
            }
         }

         if (I <= N) {
            for (J = KL + I + 1; J <= M; J++) { // 50
               A( J, I ) = ZERO;
            } // 50
         }

         if (I <= M) {
            for (J = KU + I + 1; J <= N; J++) { // 60
               A( I, J ) = ZERO;
            } // 60
         }
      } // 70
      return;
      }
