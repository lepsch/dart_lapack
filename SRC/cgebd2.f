      SUBROUTINE CGEBD2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * )
      COMPLEX            A( LDA, * ), TAUP( * ), TAUQ( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ZERO, ONE
      const              ZERO = ( 0.0E+0, 0.0E+0 ), ONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I;
      COMPLEX            ALPHA
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACGV, CLARF, CLARFG, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -4
      }
      if ( INFO.LT.0 ) {
         xerbla('CGEBD2', -INFO );
         RETURN
      }

      if ( M.GE.N ) {

         // Reduce to upper bidiagonal form

         DO 10 I = 1, N

            // Generate elementary reflector H(i) to annihilate A(i+1:m,i)

            ALPHA = A( I, I )
            clarfg(M-I+1, ALPHA, A( MIN( I+1, M ), I ), 1, TAUQ( I ) );
            D( I ) = REAL( ALPHA )
            A( I, I ) = ONE

            // Apply H(i)**H to A(i:m,i+1:n) from the left

            IF( I.LT.N ) CALL CLARF( 'Left', M-I+1, N-I, A( I, I ), 1, CONJG( TAUQ( I ) ), A( I, I+1 ), LDA, WORK )
            A( I, I ) = D( I )

            if ( I.LT.N ) {

               // Generate elementary reflector G(i) to annihilate
               // A(i,i+2:n)

               clacgv(N-I, A( I, I+1 ), LDA );
               ALPHA = A( I, I+1 )
               clarfg(N-I, ALPHA, A( I, MIN( I+2, N ) ), LDA, TAUP( I ) );
               E( I ) = REAL( ALPHA )
               A( I, I+1 ) = ONE

               // Apply G(i) to A(i+1:m,i+1:n) from the right

               clarf('Right', M-I, N-I, A( I, I+1 ), LDA, TAUP( I ), A( I+1, I+1 ), LDA, WORK );
               clacgv(N-I, A( I, I+1 ), LDA );
               A( I, I+1 ) = E( I )
            } else {
               TAUP( I ) = ZERO
            }
   10    CONTINUE
      } else {

         // Reduce to lower bidiagonal form

         DO 20 I = 1, M

            // Generate elementary reflector G(i) to annihilate A(i,i+1:n)

            clacgv(N-I+1, A( I, I ), LDA );
            ALPHA = A( I, I )
            clarfg(N-I+1, ALPHA, A( I, MIN( I+1, N ) ), LDA, TAUP( I ) );
            D( I ) = REAL( ALPHA )
            A( I, I ) = ONE

            // Apply G(i) to A(i+1:m,i:n) from the right

            IF( I.LT.M ) CALL CLARF( 'Right', M-I, N-I+1, A( I, I ), LDA, TAUP( I ), A( I+1, I ), LDA, WORK )
            clacgv(N-I+1, A( I, I ), LDA );
            A( I, I ) = D( I )

            if ( I.LT.M ) {

               // Generate elementary reflector H(i) to annihilate
               // A(i+2:m,i)

               ALPHA = A( I+1, I )
               clarfg(M-I, ALPHA, A( MIN( I+2, M ), I ), 1, TAUQ( I ) );
               E( I ) = REAL( ALPHA )
               A( I+1, I ) = ONE

               // Apply H(i)**H to A(i+1:m,i+1:n) from the left

               clarf('Left', M-I, N-I, A( I+1, I ), 1, CONJG( TAUQ( I ) ), A( I+1, I+1 ), LDA, WORK );
               A( I+1, I ) = E( I )
            } else {
               TAUQ( I ) = ZERO
            }
   20    CONTINUE
      }
      RETURN

      // End of CGEBD2

      }
