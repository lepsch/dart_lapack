      SUBROUTINE DLABRD( M, N, NB, A, LDA, D, E, TAUQ, TAUP, X, LDX, Y, LDY )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDX, LDY, M, N, NB;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), D( * ), E( * ), TAUP( * ), TAUQ( * ), X( LDX, * ), Y( LDY, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMV, DLARFG, DSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      IF( M.LE.0 .OR. N.LE.0 ) RETURN

      if ( M.GE.N ) {

         // Reduce to upper bidiagonal form

         for (I = 1; I <= NB; I++) { // 10

            // Update A(i:m,i)

            dgemv('No transpose', M-I+1, I-1, -ONE, A( I, 1 ), LDA, Y( I, 1 ), LDY, ONE, A( I, I ), 1 )             CALL DGEMV( 'No transpose', M-I+1, I-1, -ONE, X( I, 1 ), LDX, A( 1, I ), 1, ONE, A( I, I ), 1 );

            // Generate reflection Q(i) to annihilate A(i+1:m,i)

            dlarfg(M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1, TAUQ( I ) );
            D( I ) = A( I, I )
            if ( I.LT.N ) {
               A( I, I ) = ONE

               // Compute Y(i+1:n,i)

               dgemv('Transpose', M-I+1, N-I, ONE, A( I, I+1 ), LDA, A( I, I ), 1, ZERO, Y( I+1, I ), 1 )                CALL DGEMV( 'Transpose', M-I+1, I-1, ONE, A( I, 1 ), LDA, A( I, I ), 1, ZERO, Y( 1, I ), 1 )                CALL DGEMV( 'No transpose', N-I, I-1, -ONE, Y( I+1, 1 ), LDY, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 )                CALL DGEMV( 'Transpose', M-I+1, I-1, ONE, X( I, 1 ), LDX, A( I, I ), 1, ZERO, Y( 1, I ), 1 )                CALL DGEMV( 'Transpose', I-1, N-I, -ONE, A( 1, I+1 ), LDA, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 );
               dscal(N-I, TAUQ( I ), Y( I+1, I ), 1 );

               // Update A(i,i+1:n)

               dgemv('No transpose', N-I, I, -ONE, Y( I+1, 1 ), LDY, A( I, 1 ), LDA, ONE, A( I, I+1 ), LDA )                CALL DGEMV( 'Transpose', I-1, N-I, -ONE, A( 1, I+1 ), LDA, X( I, 1 ), LDX, ONE, A( I, I+1 ), LDA );

               // Generate reflection P(i) to annihilate A(i,i+2:n)

               dlarfg(N-I, A( I, I+1 ), A( I, MIN( I+2, N ) ), LDA, TAUP( I ) );
               E( I ) = A( I, I+1 )
               A( I, I+1 ) = ONE

               // Compute X(i+1:m,i)

               dgemv('No transpose', M-I, N-I, ONE, A( I+1, I+1 ), LDA, A( I, I+1 ), LDA, ZERO, X( I+1, I ), 1 )                CALL DGEMV( 'Transpose', N-I, I, ONE, Y( I+1, 1 ), LDY, A( I, I+1 ), LDA, ZERO, X( 1, I ), 1 )                CALL DGEMV( 'No transpose', M-I, I, -ONE, A( I+1, 1 ), LDA, X( 1, I ), 1, ONE, X( I+1, I ), 1 )                CALL DGEMV( 'No transpose', I-1, N-I, ONE, A( 1, I+1 ), LDA, A( I, I+1 ), LDA, ZERO, X( 1, I ), 1 )                CALL DGEMV( 'No transpose', M-I, I-1, -ONE, X( I+1, 1 ), LDX, X( 1, I ), 1, ONE, X( I+1, I ), 1 );
               dscal(M-I, TAUP( I ), X( I+1, I ), 1 );
            }
         } // 10
      } else {

         // Reduce to lower bidiagonal form

         for (I = 1; I <= NB; I++) { // 20

            // Update A(i,i:n)

            dgemv('No transpose', N-I+1, I-1, -ONE, Y( I, 1 ), LDY, A( I, 1 ), LDA, ONE, A( I, I ), LDA )             CALL DGEMV( 'Transpose', I-1, N-I+1, -ONE, A( 1, I ), LDA, X( I, 1 ), LDX, ONE, A( I, I ), LDA );

            // Generate reflection P(i) to annihilate A(i,i+1:n)

            dlarfg(N-I+1, A( I, I ), A( I, MIN( I+1, N ) ), LDA, TAUP( I ) );
            D( I ) = A( I, I )
            if ( I.LT.M ) {
               A( I, I ) = ONE

               // Compute X(i+1:m,i)

               dgemv('No transpose', M-I, N-I+1, ONE, A( I+1, I ), LDA, A( I, I ), LDA, ZERO, X( I+1, I ), 1 )                CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, Y( I, 1 ), LDY, A( I, I ), LDA, ZERO, X( 1, I ), 1 )                CALL DGEMV( 'No transpose', M-I, I-1, -ONE, A( I+1, 1 ), LDA, X( 1, I ), 1, ONE, X( I+1, I ), 1 )                CALL DGEMV( 'No transpose', I-1, N-I+1, ONE, A( 1, I ), LDA, A( I, I ), LDA, ZERO, X( 1, I ), 1 )                CALL DGEMV( 'No transpose', M-I, I-1, -ONE, X( I+1, 1 ), LDX, X( 1, I ), 1, ONE, X( I+1, I ), 1 );
               dscal(M-I, TAUP( I ), X( I+1, I ), 1 );

               // Update A(i+1:m,i)

               dgemv('No transpose', M-I, I-1, -ONE, A( I+1, 1 ), LDA, Y( I, 1 ), LDY, ONE, A( I+1, I ), 1 )                CALL DGEMV( 'No transpose', M-I, I, -ONE, X( I+1, 1 ), LDX, A( 1, I ), 1, ONE, A( I+1, I ), 1 );

               // Generate reflection Q(i) to annihilate A(i+2:m,i)

               dlarfg(M-I, A( I+1, I ), A( MIN( I+2, M ), I ), 1, TAUQ( I ) );
               E( I ) = A( I+1, I )
               A( I+1, I ) = ONE

               // Compute Y(i+1:n,i)

               dgemv('Transpose', M-I, N-I, ONE, A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, Y( I+1, I ), 1 )                CALL DGEMV( 'Transpose', M-I, I-1, ONE, A( I+1, 1 ), LDA, A( I+1, I ), 1, ZERO, Y( 1, I ), 1 )                CALL DGEMV( 'No transpose', N-I, I-1, -ONE, Y( I+1, 1 ), LDY, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 )                CALL DGEMV( 'Transpose', M-I, I, ONE, X( I+1, 1 ), LDX, A( I+1, I ), 1, ZERO, Y( 1, I ), 1 )                CALL DGEMV( 'Transpose', I, N-I, -ONE, A( 1, I+1 ), LDA, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 );
               dscal(N-I, TAUQ( I ), Y( I+1, I ), 1 );
            }
         } // 20
      }
      RETURN

      // End of DLABRD

      }
