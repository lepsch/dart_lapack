      SUBROUTINE ZLABRD( M, N, NB, A, LDA, D, E, TAUQ, TAUP, X, LDX, Y, LDY )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDX, LDY, M, N, NB;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * );
      COMPLEX*16         A( LDA, * ), TAUP( * ), TAUQ( * ), X( LDX, * ), Y( LDY, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ZERO, ONE
      const              ZERO = ( 0.0D+0, 0.0D+0 ), ONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I;
      COMPLEX*16         ALPHA
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMV, ZLACGV, ZLARFG, ZSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if (M.LE.0 || N.LE.0) RETURN;

      if ( M.GE.N ) {

         // Reduce to upper bidiagonal form

         for (I = 1; I <= NB; I++) { // 10

            // Update A(i:m,i)

            zlacgv(I-1, Y( I, 1 ), LDY );
            zgemv('No transpose', M-I+1, I-1, -ONE, A( I, 1 ), LDA, Y( I, 1 ), LDY, ONE, A( I, I ), 1 );
            zlacgv(I-1, Y( I, 1 ), LDY );
            zgemv('No transpose', M-I+1, I-1, -ONE, X( I, 1 ), LDX, A( 1, I ), 1, ONE, A( I, I ), 1 );

            // Generate reflection Q(i) to annihilate A(i+1:m,i)

            ALPHA = A( I, I )
            zlarfg(M-I+1, ALPHA, A( MIN( I+1, M ), I ), 1, TAUQ( I ) );
            D( I ) = DBLE( ALPHA )
            if ( I.LT.N ) {
               A( I, I ) = ONE

               // Compute Y(i+1:n,i)

               zgemv('Conjugate transpose', M-I+1, N-I, ONE, A( I, I+1 ), LDA, A( I, I ), 1, ZERO, Y( I+1, I ), 1 );
               zgemv('Conjugate transpose', M-I+1, I-1, ONE, A( I, 1 ), LDA, A( I, I ), 1, ZERO, Y( 1, I ), 1 );
               zgemv('No transpose', N-I, I-1, -ONE, Y( I+1, 1 ), LDY, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 );
               zgemv('Conjugate transpose', M-I+1, I-1, ONE, X( I, 1 ), LDX, A( I, I ), 1, ZERO, Y( 1, I ), 1 );
               zgemv('Conjugate transpose', I-1, N-I, -ONE, A( 1, I+1 ), LDA, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 );
               zscal(N-I, TAUQ( I ), Y( I+1, I ), 1 );

               // Update A(i,i+1:n)

               zlacgv(N-I, A( I, I+1 ), LDA );
               zlacgv(I, A( I, 1 ), LDA );
               zgemv('No transpose', N-I, I, -ONE, Y( I+1, 1 ), LDY, A( I, 1 ), LDA, ONE, A( I, I+1 ), LDA );
               zlacgv(I, A( I, 1 ), LDA );
               zlacgv(I-1, X( I, 1 ), LDX );
               zgemv('Conjugate transpose', I-1, N-I, -ONE, A( 1, I+1 ), LDA, X( I, 1 ), LDX, ONE, A( I, I+1 ), LDA );
               zlacgv(I-1, X( I, 1 ), LDX );

               // Generate reflection P(i) to annihilate A(i,i+2:n)

               ALPHA = A( I, I+1 )
               zlarfg(N-I, ALPHA, A( I, MIN( I+2, N ) ), LDA, TAUP( I ) );
               E( I ) = DBLE( ALPHA )
               A( I, I+1 ) = ONE

               // Compute X(i+1:m,i)

               zgemv('No transpose', M-I, N-I, ONE, A( I+1, I+1 ), LDA, A( I, I+1 ), LDA, ZERO, X( I+1, I ), 1 );
               zgemv('Conjugate transpose', N-I, I, ONE, Y( I+1, 1 ), LDY, A( I, I+1 ), LDA, ZERO, X( 1, I ), 1 );
               zgemv('No transpose', M-I, I, -ONE, A( I+1, 1 ), LDA, X( 1, I ), 1, ONE, X( I+1, I ), 1 );
               zgemv('No transpose', I-1, N-I, ONE, A( 1, I+1 ), LDA, A( I, I+1 ), LDA, ZERO, X( 1, I ), 1 );
               zgemv('No transpose', M-I, I-1, -ONE, X( I+1, 1 ), LDX, X( 1, I ), 1, ONE, X( I+1, I ), 1 );
               zscal(M-I, TAUP( I ), X( I+1, I ), 1 );
               zlacgv(N-I, A( I, I+1 ), LDA );
            }
         } // 10
      } else {

         // Reduce to lower bidiagonal form

         for (I = 1; I <= NB; I++) { // 20

            // Update A(i,i:n)

            zlacgv(N-I+1, A( I, I ), LDA );
            zlacgv(I-1, A( I, 1 ), LDA );
            zgemv('No transpose', N-I+1, I-1, -ONE, Y( I, 1 ), LDY, A( I, 1 ), LDA, ONE, A( I, I ), LDA );
            zlacgv(I-1, A( I, 1 ), LDA );
            zlacgv(I-1, X( I, 1 ), LDX );
            zgemv('Conjugate transpose', I-1, N-I+1, -ONE, A( 1, I ), LDA, X( I, 1 ), LDX, ONE, A( I, I ), LDA );
            zlacgv(I-1, X( I, 1 ), LDX );

            // Generate reflection P(i) to annihilate A(i,i+1:n)

            ALPHA = A( I, I )
            zlarfg(N-I+1, ALPHA, A( I, MIN( I+1, N ) ), LDA, TAUP( I ) );
            D( I ) = DBLE( ALPHA )
            if ( I.LT.M ) {
               A( I, I ) = ONE

               // Compute X(i+1:m,i)

               zgemv('No transpose', M-I, N-I+1, ONE, A( I+1, I ), LDA, A( I, I ), LDA, ZERO, X( I+1, I ), 1 );
               zgemv('Conjugate transpose', N-I+1, I-1, ONE, Y( I, 1 ), LDY, A( I, I ), LDA, ZERO, X( 1, I ), 1 );
               zgemv('No transpose', M-I, I-1, -ONE, A( I+1, 1 ), LDA, X( 1, I ), 1, ONE, X( I+1, I ), 1 );
               zgemv('No transpose', I-1, N-I+1, ONE, A( 1, I ), LDA, A( I, I ), LDA, ZERO, X( 1, I ), 1 );
               zgemv('No transpose', M-I, I-1, -ONE, X( I+1, 1 ), LDX, X( 1, I ), 1, ONE, X( I+1, I ), 1 );
               zscal(M-I, TAUP( I ), X( I+1, I ), 1 );
               zlacgv(N-I+1, A( I, I ), LDA );

               // Update A(i+1:m,i)

               zlacgv(I-1, Y( I, 1 ), LDY );
               zgemv('No transpose', M-I, I-1, -ONE, A( I+1, 1 ), LDA, Y( I, 1 ), LDY, ONE, A( I+1, I ), 1 );
               zlacgv(I-1, Y( I, 1 ), LDY );
               zgemv('No transpose', M-I, I, -ONE, X( I+1, 1 ), LDX, A( 1, I ), 1, ONE, A( I+1, I ), 1 );

               // Generate reflection Q(i) to annihilate A(i+2:m,i)

               ALPHA = A( I+1, I )
               zlarfg(M-I, ALPHA, A( MIN( I+2, M ), I ), 1, TAUQ( I ) );
               E( I ) = DBLE( ALPHA )
               A( I+1, I ) = ONE

               // Compute Y(i+1:n,i)

               zgemv('Conjugate transpose', M-I, N-I, ONE, A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, Y( I+1, I ), 1 );
               zgemv('Conjugate transpose', M-I, I-1, ONE, A( I+1, 1 ), LDA, A( I+1, I ), 1, ZERO, Y( 1, I ), 1 );
               zgemv('No transpose', N-I, I-1, -ONE, Y( I+1, 1 ), LDY, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 );
               zgemv('Conjugate transpose', M-I, I, ONE, X( I+1, 1 ), LDX, A( I+1, I ), 1, ZERO, Y( 1, I ), 1 );
               zgemv('Conjugate transpose', I, N-I, -ONE, A( 1, I+1 ), LDA, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 );
               zscal(N-I, TAUQ( I ), Y( I+1, I ), 1 );
            } else {
               zlacgv(N-I+1, A( I, I ), LDA );
            }
         } // 20
      }
      RETURN

      // End of ZLABRD

      }
