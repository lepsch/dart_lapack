      SUBROUTINE SLABRD( M, N, NB, A, LDA, D, E, TAUQ, TAUP, X, LDX, Y, LDY )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDX, LDY, M, N, NB;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), D( * ), E( * ), TAUP( * ), TAUQ( * ), X( LDX, * ), Y( LDY, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMV, SLARFG, SSCAL
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

            sgemv('No transpose', M-I+1, I-1, -ONE, A( I, 1 ), LDA, Y( I, 1 ), LDY, ONE, A( I, I ), 1 )             CALL SGEMV( 'No transpose', M-I+1, I-1, -ONE, X( I, 1 ), LDX, A( 1, I ), 1, ONE, A( I, I ), 1 );

            // Generate reflection Q(i) to annihilate A(i+1:m,i)

            slarfg(M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1, TAUQ( I ) );
            D( I ) = A( I, I )
            if ( I.LT.N ) {
               A( I, I ) = ONE

               // Compute Y(i+1:n,i)

               sgemv('Transpose', M-I+1, N-I, ONE, A( I, I+1 ), LDA, A( I, I ), 1, ZERO, Y( I+1, I ), 1 )                CALL SGEMV( 'Transpose', M-I+1, I-1, ONE, A( I, 1 ), LDA, A( I, I ), 1, ZERO, Y( 1, I ), 1 )                CALL SGEMV( 'No transpose', N-I, I-1, -ONE, Y( I+1, 1 ), LDY, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 )                CALL SGEMV( 'Transpose', M-I+1, I-1, ONE, X( I, 1 ), LDX, A( I, I ), 1, ZERO, Y( 1, I ), 1 )                CALL SGEMV( 'Transpose', I-1, N-I, -ONE, A( 1, I+1 ), LDA, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 );
               sscal(N-I, TAUQ( I ), Y( I+1, I ), 1 );

               // Update A(i,i+1:n)

               sgemv('No transpose', N-I, I, -ONE, Y( I+1, 1 ), LDY, A( I, 1 ), LDA, ONE, A( I, I+1 ), LDA )                CALL SGEMV( 'Transpose', I-1, N-I, -ONE, A( 1, I+1 ), LDA, X( I, 1 ), LDX, ONE, A( I, I+1 ), LDA );

               // Generate reflection P(i) to annihilate A(i,i+2:n)

               slarfg(N-I, A( I, I+1 ), A( I, MIN( I+2, N ) ), LDA, TAUP( I ) );
               E( I ) = A( I, I+1 )
               A( I, I+1 ) = ONE

               // Compute X(i+1:m,i)

               sgemv('No transpose', M-I, N-I, ONE, A( I+1, I+1 ), LDA, A( I, I+1 ), LDA, ZERO, X( I+1, I ), 1 )                CALL SGEMV( 'Transpose', N-I, I, ONE, Y( I+1, 1 ), LDY, A( I, I+1 ), LDA, ZERO, X( 1, I ), 1 )                CALL SGEMV( 'No transpose', M-I, I, -ONE, A( I+1, 1 ), LDA, X( 1, I ), 1, ONE, X( I+1, I ), 1 )                CALL SGEMV( 'No transpose', I-1, N-I, ONE, A( 1, I+1 ), LDA, A( I, I+1 ), LDA, ZERO, X( 1, I ), 1 )                CALL SGEMV( 'No transpose', M-I, I-1, -ONE, X( I+1, 1 ), LDX, X( 1, I ), 1, ONE, X( I+1, I ), 1 );
               sscal(M-I, TAUP( I ), X( I+1, I ), 1 );
            }
         } // 10
      } else {

         // Reduce to lower bidiagonal form

         for (I = 1; I <= NB; I++) { // 20

            // Update A(i,i:n)

            sgemv('No transpose', N-I+1, I-1, -ONE, Y( I, 1 ), LDY, A( I, 1 ), LDA, ONE, A( I, I ), LDA )             CALL SGEMV( 'Transpose', I-1, N-I+1, -ONE, A( 1, I ), LDA, X( I, 1 ), LDX, ONE, A( I, I ), LDA );

            // Generate reflection P(i) to annihilate A(i,i+1:n)

            slarfg(N-I+1, A( I, I ), A( I, MIN( I+1, N ) ), LDA, TAUP( I ) );
            D( I ) = A( I, I )
            if ( I.LT.M ) {
               A( I, I ) = ONE

               // Compute X(i+1:m,i)

               sgemv('No transpose', M-I, N-I+1, ONE, A( I+1, I ), LDA, A( I, I ), LDA, ZERO, X( I+1, I ), 1 )                CALL SGEMV( 'Transpose', N-I+1, I-1, ONE, Y( I, 1 ), LDY, A( I, I ), LDA, ZERO, X( 1, I ), 1 )                CALL SGEMV( 'No transpose', M-I, I-1, -ONE, A( I+1, 1 ), LDA, X( 1, I ), 1, ONE, X( I+1, I ), 1 )                CALL SGEMV( 'No transpose', I-1, N-I+1, ONE, A( 1, I ), LDA, A( I, I ), LDA, ZERO, X( 1, I ), 1 )                CALL SGEMV( 'No transpose', M-I, I-1, -ONE, X( I+1, 1 ), LDX, X( 1, I ), 1, ONE, X( I+1, I ), 1 );
               sscal(M-I, TAUP( I ), X( I+1, I ), 1 );

               // Update A(i+1:m,i)

               sgemv('No transpose', M-I, I-1, -ONE, A( I+1, 1 ), LDA, Y( I, 1 ), LDY, ONE, A( I+1, I ), 1 )                CALL SGEMV( 'No transpose', M-I, I, -ONE, X( I+1, 1 ), LDX, A( 1, I ), 1, ONE, A( I+1, I ), 1 );

               // Generate reflection Q(i) to annihilate A(i+2:m,i)

               slarfg(M-I, A( I+1, I ), A( MIN( I+2, M ), I ), 1, TAUQ( I ) );
               E( I ) = A( I+1, I )
               A( I+1, I ) = ONE

               // Compute Y(i+1:n,i)

               sgemv('Transpose', M-I, N-I, ONE, A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, Y( I+1, I ), 1 )                CALL SGEMV( 'Transpose', M-I, I-1, ONE, A( I+1, 1 ), LDA, A( I+1, I ), 1, ZERO, Y( 1, I ), 1 )                CALL SGEMV( 'No transpose', N-I, I-1, -ONE, Y( I+1, 1 ), LDY, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 )                CALL SGEMV( 'Transpose', M-I, I, ONE, X( I+1, 1 ), LDX, A( I+1, I ), 1, ZERO, Y( 1, I ), 1 )                CALL SGEMV( 'Transpose', I, N-I, -ONE, A( 1, I+1 ), LDA, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 );
               sscal(N-I, TAUQ( I ), Y( I+1, I ), 1 );
            }
         } // 20
      }
      RETURN

      // End of SLABRD

      }
