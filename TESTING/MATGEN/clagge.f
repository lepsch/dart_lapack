      SUBROUTINE CLAGGE( M, N, KL, KU, D, A, LDA, ISEED, WORK, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, KL, KU, LDA, M, N;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      REAL               D( * )
      COMPLEX            A( LDA, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ZERO, ONE
      const              ZERO = ( 0.0E+0, 0.0E+0 ), ONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               WN
      COMPLEX            TAU, WA, WB
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV, CGERC, CLACGV, CLARNV, CSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL
      // ..
      // .. External Functions ..
      REAL               SCNRM2
      // EXTERNAL SCNRM2
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      if ( M < 0 ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( KL < 0 || KL.GT.M-1 ) {
         INFO = -3
      } else if ( KU < 0 || KU.GT.N-1 ) {
         INFO = -4
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -7
      }
      if ( INFO < 0 ) {
         xerbla('CLAGGE', -INFO );
         RETURN
      }

      // initialize A to diagonal matrix

      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= M; I++) { // 10
            A( I, J ) = ZERO
         } // 10
      } // 20
      DO 30 I = 1, MIN( M, N )
         A( I, I ) = D( I )
      } // 30

      // Quick exit if the user wants a diagonal matrix

      IF(( KL == 0 ) && ( KU == 0)) RETURN

      // pre- and post-multiply A by random unitary matrices

      DO 40 I = MIN( M, N ), 1, -1
         if ( I < M ) {

            // generate random reflection

            clarnv(3, ISEED, M-I+1, WORK );
            WN = SCNRM2( M-I+1, WORK, 1 )
            WA = ( WN / ABS( WORK( 1 ) ) )*WORK( 1 )
            if ( WN == ZERO ) {
               TAU = ZERO
            } else {
               WB = WORK( 1 ) + WA
               cscal(M-I, ONE / WB, WORK( 2 ), 1 );
               WORK( 1 ) = ONE
               TAU = REAL( WB / WA )
            }

            // multiply A(i:m,i:n) by random reflection from the left

            cgemv('Conjugate transpose', M-I+1, N-I+1, ONE, A( I, I ), LDA, WORK, 1, ZERO, WORK( M+1 ), 1 );
            cgerc(M-I+1, N-I+1, -TAU, WORK, 1, WORK( M+1 ), 1, A( I, I ), LDA );
         }
         if ( I < N ) {

            // generate random reflection

            clarnv(3, ISEED, N-I+1, WORK );
            WN = SCNRM2( N-I+1, WORK, 1 )
            WA = ( WN / ABS( WORK( 1 ) ) )*WORK( 1 )
            if ( WN == ZERO ) {
               TAU = ZERO
            } else {
               WB = WORK( 1 ) + WA
               cscal(N-I, ONE / WB, WORK( 2 ), 1 );
               WORK( 1 ) = ONE
               TAU = REAL( WB / WA )
            }

            // multiply A(i:m,i:n) by random reflection from the right

            cgemv('No transpose', M-I+1, N-I+1, ONE, A( I, I ), LDA, WORK, 1, ZERO, WORK( N+1 ), 1 );
            cgerc(M-I+1, N-I+1, -TAU, WORK( N+1 ), 1, WORK, 1, A( I, I ), LDA );
         }
      } // 40

      // Reduce number of subdiagonals to KL and number of superdiagonals
      // to KU

      DO 70 I = 1, MAX( M-1-KL, N-1-KU )
         if ( KL.LE.KU ) {

            // annihilate subdiagonal elements first (necessary if KL = 0)

            if ( I.LE.MIN( M-1-KL, N ) ) {

               // generate reflection to annihilate A(kl+i+1:m,i)

               WN = SCNRM2( M-KL-I+1, A( KL+I, I ), 1 )
               WA = ( WN / ABS( A( KL+I, I ) ) )*A( KL+I, I )
               if ( WN == ZERO ) {
                  TAU = ZERO
               } else {
                  WB = A( KL+I, I ) + WA
                  cscal(M-KL-I, ONE / WB, A( KL+I+1, I ), 1 );
                  A( KL+I, I ) = ONE
                  TAU = REAL( WB / WA )
               }

               // apply reflection to A(kl+i:m,i+1:n) from the left

               cgemv('Conjugate transpose', M-KL-I+1, N-I, ONE, A( KL+I, I+1 ), LDA, A( KL+I, I ), 1, ZERO, WORK, 1 );
               cgerc(M-KL-I+1, N-I, -TAU, A( KL+I, I ), 1, WORK, 1, A( KL+I, I+1 ), LDA );
               A( KL+I, I ) = -WA
            }

            if ( I.LE.MIN( N-1-KU, M ) ) {

               // generate reflection to annihilate A(i,ku+i+1:n)

               WN = SCNRM2( N-KU-I+1, A( I, KU+I ), LDA )
               WA = ( WN / ABS( A( I, KU+I ) ) )*A( I, KU+I )
               if ( WN == ZERO ) {
                  TAU = ZERO
               } else {
                  WB = A( I, KU+I ) + WA
                  cscal(N-KU-I, ONE / WB, A( I, KU+I+1 ), LDA );
                  A( I, KU+I ) = ONE
                  TAU = REAL( WB / WA )
               }

               // apply reflection to A(i+1:m,ku+i:n) from the right

               clacgv(N-KU-I+1, A( I, KU+I ), LDA );
               cgemv('No transpose', M-I, N-KU-I+1, ONE, A( I+1, KU+I ), LDA, A( I, KU+I ), LDA, ZERO, WORK, 1 );
               cgerc(M-I, N-KU-I+1, -TAU, WORK, 1, A( I, KU+I ), LDA, A( I+1, KU+I ), LDA );
               A( I, KU+I ) = -WA
            }
         } else {

            // annihilate superdiagonal elements first (necessary if
            // KU = 0)

            if ( I.LE.MIN( N-1-KU, M ) ) {

               // generate reflection to annihilate A(i,ku+i+1:n)

               WN = SCNRM2( N-KU-I+1, A( I, KU+I ), LDA )
               WA = ( WN / ABS( A( I, KU+I ) ) )*A( I, KU+I )
               if ( WN == ZERO ) {
                  TAU = ZERO
               } else {
                  WB = A( I, KU+I ) + WA
                  cscal(N-KU-I, ONE / WB, A( I, KU+I+1 ), LDA );
                  A( I, KU+I ) = ONE
                  TAU = REAL( WB / WA )
               }

               // apply reflection to A(i+1:m,ku+i:n) from the right

               clacgv(N-KU-I+1, A( I, KU+I ), LDA );
               cgemv('No transpose', M-I, N-KU-I+1, ONE, A( I+1, KU+I ), LDA, A( I, KU+I ), LDA, ZERO, WORK, 1 );
               cgerc(M-I, N-KU-I+1, -TAU, WORK, 1, A( I, KU+I ), LDA, A( I+1, KU+I ), LDA );
               A( I, KU+I ) = -WA
            }

            if ( I.LE.MIN( M-1-KL, N ) ) {

               // generate reflection to annihilate A(kl+i+1:m,i)

               WN = SCNRM2( M-KL-I+1, A( KL+I, I ), 1 )
               WA = ( WN / ABS( A( KL+I, I ) ) )*A( KL+I, I )
               if ( WN == ZERO ) {
                  TAU = ZERO
               } else {
                  WB = A( KL+I, I ) + WA
                  cscal(M-KL-I, ONE / WB, A( KL+I+1, I ), 1 );
                  A( KL+I, I ) = ONE
                  TAU = REAL( WB / WA )
               }

               // apply reflection to A(kl+i:m,i+1:n) from the left

               cgemv('Conjugate transpose', M-KL-I+1, N-I, ONE, A( KL+I, I+1 ), LDA, A( KL+I, I ), 1, ZERO, WORK, 1 );
               cgerc(M-KL-I+1, N-I, -TAU, A( KL+I, I ), 1, WORK, 1, A( KL+I, I+1 ), LDA );
               A( KL+I, I ) = -WA
            }
         }

         if (I .LE. N) {
            for (J = KL + I + 1; J <= M; J++) { // 50
               A( J, I ) = ZERO
            } // 50
         }

         if (I .LE. M) {
            for (J = KU + I + 1; J <= N; J++) { // 60
               A( I, J ) = ZERO
            } // 60
         }
      } // 70
      RETURN

      // End of CLAGGE

      }
