      SUBROUTINE DLAGGE( M, N, KL, KU, D, A, LDA, ISEED, WORK, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, KL, KU, LDA, M, N;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      double             A( LDA, * ), D( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
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
      double             DNRM2;
      // EXTERNAL DNRM2
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( KL.LT.0 .OR. KL.GT.M-1 ) {
         INFO = -3
      } else if ( KU.LT.0 .OR. KU.GT.N-1 ) {
         INFO = -4
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -7
      }
      if ( INFO.LT.0 ) {
         CALL XERBLA( 'DLAGGE', -INFO )
         RETURN
      }

      // initialize A to diagonal matrix

      DO 20 J = 1, N
         DO 10 I = 1, M
            A( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
      DO 30 I = 1, MIN( M, N )
         A( I, I ) = D( I )
   30 CONTINUE

      // Quick exit if the user wants a diagonal matrix

      IF(( KL .EQ. 0 ).AND.( KU .EQ. 0)) RETURN

      // pre- and post-multiply A by random orthogonal matrices

      DO 40 I = MIN( M, N ), 1, -1
         if ( I.LT.M ) {

            // generate random reflection

            CALL DLARNV( 3, ISEED, M-I+1, WORK )
            WN = DNRM2( M-I+1, WORK, 1 )
            WA = SIGN( WN, WORK( 1 ) )
            if ( WN.EQ.ZERO ) {
               TAU = ZERO
            } else {
               WB = WORK( 1 ) + WA
               CALL DSCAL( M-I, ONE / WB, WORK( 2 ), 1 )
               WORK( 1 ) = ONE
               TAU = WB / WA
            }

            // multiply A(i:m,i:n) by random reflection from the left

            CALL DGEMV( 'Transpose', M-I+1, N-I+1, ONE, A( I, I ), LDA, WORK, 1, ZERO, WORK( M+1 ), 1 )             CALL DGER( M-I+1, N-I+1, -TAU, WORK, 1, WORK( M+1 ), 1, A( I, I ), LDA )
         }
         if ( I.LT.N ) {

            // generate random reflection

            CALL DLARNV( 3, ISEED, N-I+1, WORK )
            WN = DNRM2( N-I+1, WORK, 1 )
            WA = SIGN( WN, WORK( 1 ) )
            if ( WN.EQ.ZERO ) {
               TAU = ZERO
            } else {
               WB = WORK( 1 ) + WA
               CALL DSCAL( N-I, ONE / WB, WORK( 2 ), 1 )
               WORK( 1 ) = ONE
               TAU = WB / WA
            }

            // multiply A(i:m,i:n) by random reflection from the right

            CALL DGEMV( 'No transpose', M-I+1, N-I+1, ONE, A( I, I ), LDA, WORK, 1, ZERO, WORK( N+1 ), 1 )             CALL DGER( M-I+1, N-I+1, -TAU, WORK( N+1 ), 1, WORK, 1, A( I, I ), LDA )
         }
   40 CONTINUE

      // Reduce number of subdiagonals to KL and number of superdiagonals
     t // o KU

      DO 70 I = 1, MAX( M-1-KL, N-1-KU )
         if ( KL.LE.KU ) {

            // annihilate subdiagonal elements first (necessary if KL = 0)

            if ( I.LE.MIN( M-1-KL, N ) ) {

               // generate reflection to annihilate A(kl+i+1:m,i)

               WN = DNRM2( M-KL-I+1, A( KL+I, I ), 1 )
               WA = SIGN( WN, A( KL+I, I ) )
               if ( WN.EQ.ZERO ) {
                  TAU = ZERO
               } else {
                  WB = A( KL+I, I ) + WA
                  CALL DSCAL( M-KL-I, ONE / WB, A( KL+I+1, I ), 1 )
                  A( KL+I, I ) = ONE
                  TAU = WB / WA
               }

               // apply reflection to A(kl+i:m,i+1:n) from the left

               CALL DGEMV( 'Transpose', M-KL-I+1, N-I, ONE, A( KL+I, I+1 ), LDA, A( KL+I, I ), 1, ZERO, WORK, 1 )
               CALL DGER( M-KL-I+1, N-I, -TAU, A( KL+I, I ), 1, WORK, 1, A( KL+I, I+1 ), LDA )
               A( KL+I, I ) = -WA
            }

            if ( I.LE.MIN( N-1-KU, M ) ) {

               // generate reflection to annihilate A(i,ku+i+1:n)

               WN = DNRM2( N-KU-I+1, A( I, KU+I ), LDA )
               WA = SIGN( WN, A( I, KU+I ) )
               if ( WN.EQ.ZERO ) {
                  TAU = ZERO
               } else {
                  WB = A( I, KU+I ) + WA
                  CALL DSCAL( N-KU-I, ONE / WB, A( I, KU+I+1 ), LDA )
                  A( I, KU+I ) = ONE
                  TAU = WB / WA
               }

               // apply reflection to A(i+1:m,ku+i:n) from the right

               CALL DGEMV( 'No transpose', M-I, N-KU-I+1, ONE, A( I+1, KU+I ), LDA, A( I, KU+I ), LDA, ZERO, WORK, 1 )
               CALL DGER( M-I, N-KU-I+1, -TAU, WORK, 1, A( I, KU+I ), LDA, A( I+1, KU+I ), LDA )
               A( I, KU+I ) = -WA
            }
         } else {

            // annihilate superdiagonal elements first (necessary if
            // KU = 0)

            if ( I.LE.MIN( N-1-KU, M ) ) {

               // generate reflection to annihilate A(i,ku+i+1:n)

               WN = DNRM2( N-KU-I+1, A( I, KU+I ), LDA )
               WA = SIGN( WN, A( I, KU+I ) )
               if ( WN.EQ.ZERO ) {
                  TAU = ZERO
               } else {
                  WB = A( I, KU+I ) + WA
                  CALL DSCAL( N-KU-I, ONE / WB, A( I, KU+I+1 ), LDA )
                  A( I, KU+I ) = ONE
                  TAU = WB / WA
               }

               // apply reflection to A(i+1:m,ku+i:n) from the right

               CALL DGEMV( 'No transpose', M-I, N-KU-I+1, ONE, A( I+1, KU+I ), LDA, A( I, KU+I ), LDA, ZERO, WORK, 1 )
               CALL DGER( M-I, N-KU-I+1, -TAU, WORK, 1, A( I, KU+I ), LDA, A( I+1, KU+I ), LDA )
               A( I, KU+I ) = -WA
            }

            if ( I.LE.MIN( M-1-KL, N ) ) {

               // generate reflection to annihilate A(kl+i+1:m,i)

               WN = DNRM2( M-KL-I+1, A( KL+I, I ), 1 )
               WA = SIGN( WN, A( KL+I, I ) )
               if ( WN.EQ.ZERO ) {
                  TAU = ZERO
               } else {
                  WB = A( KL+I, I ) + WA
                  CALL DSCAL( M-KL-I, ONE / WB, A( KL+I+1, I ), 1 )
                  A( KL+I, I ) = ONE
                  TAU = WB / WA
               }

               // apply reflection to A(kl+i:m,i+1:n) from the left

               CALL DGEMV( 'Transpose', M-KL-I+1, N-I, ONE, A( KL+I, I+1 ), LDA, A( KL+I, I ), 1, ZERO, WORK, 1 )
               CALL DGER( M-KL-I+1, N-I, -TAU, A( KL+I, I ), 1, WORK, 1, A( KL+I, I+1 ), LDA )
               A( KL+I, I ) = -WA
            }
         }

         if (I .LE. N) {
            DO 50 J = KL + I + 1, M
               A( J, I ) = ZERO
   50       CONTINUE
         }

         if (I .LE. M) {
            DO 60 J = KU + I + 1, N
               A( I, J ) = ZERO
   60       CONTINUE
         }
   70 CONTINUE
      RETURN

      // End of DLAGGE

      }
