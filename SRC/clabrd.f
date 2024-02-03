      SUBROUTINE CLABRD( M, N, NB, A, LDA, D, E, TAUQ, TAUP, X, LDX, Y, LDY );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDX, LDY, M, N, NB;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * );
      COMPLEX            A( LDA, * ), TAUP( * ), TAUQ( * ), X( LDX, * ), Y( LDY, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I;
      COMPLEX            ALPHA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV, CLACGV, CLARFG, CSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if (M <= 0 || N <= 0) RETURN;

      if ( M >= N ) {

         // Reduce to upper bidiagonal form

         for (I = 1; I <= NB; I++) { // 10

            // Update A(i:m,i)

            clacgv(I-1, Y( I, 1 ), LDY );
            cgemv('No transpose', M-I+1, I-1, -ONE, A( I, 1 ), LDA, Y( I, 1 ), LDY, ONE, A( I, I ), 1 );
            clacgv(I-1, Y( I, 1 ), LDY );
            cgemv('No transpose', M-I+1, I-1, -ONE, X( I, 1 ), LDX, A( 1, I ), 1, ONE, A( I, I ), 1 );

            // Generate reflection Q(i) to annihilate A(i+1:m,i)

            ALPHA = A( I, I );
            clarfg(M-I+1, ALPHA, A( MIN( I+1, M ), I ), 1, TAUQ( I ) );
            D( I ) = REAL( ALPHA );
            if ( I < N ) {
               A( I, I ) = ONE;

               // Compute Y(i+1:n,i)

               cgemv('Conjugate transpose', M-I+1, N-I, ONE, A( I, I+1 ), LDA, A( I, I ), 1, ZERO, Y( I+1, I ), 1 );
               cgemv('Conjugate transpose', M-I+1, I-1, ONE, A( I, 1 ), LDA, A( I, I ), 1, ZERO, Y( 1, I ), 1 );
               cgemv('No transpose', N-I, I-1, -ONE, Y( I+1, 1 ), LDY, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 );
               cgemv('Conjugate transpose', M-I+1, I-1, ONE, X( I, 1 ), LDX, A( I, I ), 1, ZERO, Y( 1, I ), 1 );
               cgemv('Conjugate transpose', I-1, N-I, -ONE, A( 1, I+1 ), LDA, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 );
               cscal(N-I, TAUQ( I ), Y( I+1, I ), 1 );

               // Update A(i,i+1:n)

               clacgv(N-I, A( I, I+1 ), LDA );
               clacgv(I, A( I, 1 ), LDA );
               cgemv('No transpose', N-I, I, -ONE, Y( I+1, 1 ), LDY, A( I, 1 ), LDA, ONE, A( I, I+1 ), LDA );
               clacgv(I, A( I, 1 ), LDA );
               clacgv(I-1, X( I, 1 ), LDX );
               cgemv('Conjugate transpose', I-1, N-I, -ONE, A( 1, I+1 ), LDA, X( I, 1 ), LDX, ONE, A( I, I+1 ), LDA );
               clacgv(I-1, X( I, 1 ), LDX );

               // Generate reflection P(i) to annihilate A(i,i+2:n)

               ALPHA = A( I, I+1 );
               clarfg(N-I, ALPHA, A( I, MIN( I+2, N ) ), LDA, TAUP( I ) );
               E( I ) = REAL( ALPHA );
               A( I, I+1 ) = ONE;

               // Compute X(i+1:m,i)

               cgemv('No transpose', M-I, N-I, ONE, A( I+1, I+1 ), LDA, A( I, I+1 ), LDA, ZERO, X( I+1, I ), 1 );
               cgemv('Conjugate transpose', N-I, I, ONE, Y( I+1, 1 ), LDY, A( I, I+1 ), LDA, ZERO, X( 1, I ), 1 );
               cgemv('No transpose', M-I, I, -ONE, A( I+1, 1 ), LDA, X( 1, I ), 1, ONE, X( I+1, I ), 1 );
               cgemv('No transpose', I-1, N-I, ONE, A( 1, I+1 ), LDA, A( I, I+1 ), LDA, ZERO, X( 1, I ), 1 );
               cgemv('No transpose', M-I, I-1, -ONE, X( I+1, 1 ), LDX, X( 1, I ), 1, ONE, X( I+1, I ), 1 );
               cscal(M-I, TAUP( I ), X( I+1, I ), 1 );
               clacgv(N-I, A( I, I+1 ), LDA );
            }
         } // 10
      } else {

         // Reduce to lower bidiagonal form

         for (I = 1; I <= NB; I++) { // 20

            // Update A(i,i:n)

            clacgv(N-I+1, A( I, I ), LDA );
            clacgv(I-1, A( I, 1 ), LDA );
            cgemv('No transpose', N-I+1, I-1, -ONE, Y( I, 1 ), LDY, A( I, 1 ), LDA, ONE, A( I, I ), LDA );
            clacgv(I-1, A( I, 1 ), LDA );
            clacgv(I-1, X( I, 1 ), LDX );
            cgemv('Conjugate transpose', I-1, N-I+1, -ONE, A( 1, I ), LDA, X( I, 1 ), LDX, ONE, A( I, I ), LDA );
            clacgv(I-1, X( I, 1 ), LDX );

            // Generate reflection P(i) to annihilate A(i,i+1:n)

            ALPHA = A( I, I );
            clarfg(N-I+1, ALPHA, A( I, MIN( I+1, N ) ), LDA, TAUP( I ) );
            D( I ) = REAL( ALPHA );
            if ( I < M ) {
               A( I, I ) = ONE;

               // Compute X(i+1:m,i)

               cgemv('No transpose', M-I, N-I+1, ONE, A( I+1, I ), LDA, A( I, I ), LDA, ZERO, X( I+1, I ), 1 );
               cgemv('Conjugate transpose', N-I+1, I-1, ONE, Y( I, 1 ), LDY, A( I, I ), LDA, ZERO, X( 1, I ), 1 );
               cgemv('No transpose', M-I, I-1, -ONE, A( I+1, 1 ), LDA, X( 1, I ), 1, ONE, X( I+1, I ), 1 );
               cgemv('No transpose', I-1, N-I+1, ONE, A( 1, I ), LDA, A( I, I ), LDA, ZERO, X( 1, I ), 1 );
               cgemv('No transpose', M-I, I-1, -ONE, X( I+1, 1 ), LDX, X( 1, I ), 1, ONE, X( I+1, I ), 1 );
               cscal(M-I, TAUP( I ), X( I+1, I ), 1 );
               clacgv(N-I+1, A( I, I ), LDA );

               // Update A(i+1:m,i)

               clacgv(I-1, Y( I, 1 ), LDY );
               cgemv('No transpose', M-I, I-1, -ONE, A( I+1, 1 ), LDA, Y( I, 1 ), LDY, ONE, A( I+1, I ), 1 );
               clacgv(I-1, Y( I, 1 ), LDY );
               cgemv('No transpose', M-I, I, -ONE, X( I+1, 1 ), LDX, A( 1, I ), 1, ONE, A( I+1, I ), 1 );

               // Generate reflection Q(i) to annihilate A(i+2:m,i)

               ALPHA = A( I+1, I );
               clarfg(M-I, ALPHA, A( MIN( I+2, M ), I ), 1, TAUQ( I ) );
               E( I ) = REAL( ALPHA );
               A( I+1, I ) = ONE;

               // Compute Y(i+1:n,i)

               cgemv('Conjugate transpose', M-I, N-I, ONE, A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, Y( I+1, I ), 1 );
               cgemv('Conjugate transpose', M-I, I-1, ONE, A( I+1, 1 ), LDA, A( I+1, I ), 1, ZERO, Y( 1, I ), 1 );
               cgemv('No transpose', N-I, I-1, -ONE, Y( I+1, 1 ), LDY, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 );
               cgemv('Conjugate transpose', M-I, I, ONE, X( I+1, 1 ), LDX, A( I+1, I ), 1, ZERO, Y( 1, I ), 1 );
               cgemv('Conjugate transpose', I, N-I, -ONE, A( 1, I+1 ), LDA, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 );
               cscal(N-I, TAUQ( I ), Y( I+1, I ), 1 );
            } else {
               clacgv(N-I+1, A( I, I ), LDA );
            }
         } // 20
      }
      return;

      // End of CLABRD

      }
