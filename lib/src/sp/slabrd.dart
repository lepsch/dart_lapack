      void slabrd(final int M, final int N, final int NB, final Matrix<double> A, final int LDA, final int D, final int E, final int TAUQ, final int TAUP, final Matrix<double> X, final int LDX, final int Y, final int LDY,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDX, LDY, M, N, NB;
      double               A( LDA, * ), D( * ), E( * ), TAUP( * ), TAUQ( * ), X( LDX, * ), Y( LDY, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMV, SLARFG, SSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN

      // Quick return if possible

      if (M <= 0 || N <= 0) return;

      if ( M >= N ) {

         // Reduce to upper bidiagonal form

         for (I = 1; I <= NB; I++) { // 10

            // Update A(i:m,i)

            sgemv('No transpose', M-I+1, I-1, -ONE, A( I, 1 ), LDA, Y( I, 1 ), LDY, ONE, A( I, I ), 1 );
            sgemv('No transpose', M-I+1, I-1, -ONE, X( I, 1 ), LDX, A( 1, I ), 1, ONE, A( I, I ), 1 );

            // Generate reflection Q(i) to annihilate A(i+1:m,i)

            slarfg(M-I+1, A( I, I ), A( min( I+1, M ), I ), 1, TAUQ( I ) );
            D[I] = A( I, I );
            if ( I < N ) {
               A[I][I] = ONE;

               // Compute Y(i+1:n,i)

               sgemv('Transpose', M-I+1, N-I, ONE, A( I, I+1 ), LDA, A( I, I ), 1, ZERO, Y( I+1, I ), 1 );
               sgemv('Transpose', M-I+1, I-1, ONE, A( I, 1 ), LDA, A( I, I ), 1, ZERO, Y( 1, I ), 1 );
               sgemv('No transpose', N-I, I-1, -ONE, Y( I+1, 1 ), LDY, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 );
               sgemv('Transpose', M-I+1, I-1, ONE, X( I, 1 ), LDX, A( I, I ), 1, ZERO, Y( 1, I ), 1 );
               sgemv('Transpose', I-1, N-I, -ONE, A( 1, I+1 ), LDA, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 );
               sscal(N-I, TAUQ( I ), Y( I+1, I ), 1 );

               // Update A(i,i+1:n)

               sgemv('No transpose', N-I, I, -ONE, Y( I+1, 1 ), LDY, A( I, 1 ), LDA, ONE, A( I, I+1 ), LDA );
               sgemv('Transpose', I-1, N-I, -ONE, A( 1, I+1 ), LDA, X( I, 1 ), LDX, ONE, A( I, I+1 ), LDA );

               // Generate reflection P(i) to annihilate A(i,i+2:n)

               slarfg(N-I, A( I, I+1 ), A( I, min( I+2, N ) ), LDA, TAUP( I ) );
               E[I] = A( I, I+1 );
               A[I][I+1] = ONE;

               // Compute X(i+1:m,i)

               sgemv('No transpose', M-I, N-I, ONE, A( I+1, I+1 ), LDA, A( I, I+1 ), LDA, ZERO, X( I+1, I ), 1 );
               sgemv('Transpose', N-I, I, ONE, Y( I+1, 1 ), LDY, A( I, I+1 ), LDA, ZERO, X( 1, I ), 1 );
               sgemv('No transpose', M-I, I, -ONE, A( I+1, 1 ), LDA, X( 1, I ), 1, ONE, X( I+1, I ), 1 );
               sgemv('No transpose', I-1, N-I, ONE, A( 1, I+1 ), LDA, A( I, I+1 ), LDA, ZERO, X( 1, I ), 1 );
               sgemv('No transpose', M-I, I-1, -ONE, X( I+1, 1 ), LDX, X( 1, I ), 1, ONE, X( I+1, I ), 1 );
               sscal(M-I, TAUP( I ), X( I+1, I ), 1 );
            }
         } // 10
      } else {

         // Reduce to lower bidiagonal form

         for (I = 1; I <= NB; I++) { // 20

            // Update A(i,i:n)

            sgemv('No transpose', N-I+1, I-1, -ONE, Y( I, 1 ), LDY, A( I, 1 ), LDA, ONE, A( I, I ), LDA );
            sgemv('Transpose', I-1, N-I+1, -ONE, A( 1, I ), LDA, X( I, 1 ), LDX, ONE, A( I, I ), LDA );

            // Generate reflection P(i) to annihilate A(i,i+1:n)

            slarfg(N-I+1, A( I, I ), A( I, min( I+1, N ) ), LDA, TAUP( I ) );
            D[I] = A( I, I );
            if ( I < M ) {
               A[I][I] = ONE;

               // Compute X(i+1:m,i)

               sgemv('No transpose', M-I, N-I+1, ONE, A( I+1, I ), LDA, A( I, I ), LDA, ZERO, X( I+1, I ), 1 );
               sgemv('Transpose', N-I+1, I-1, ONE, Y( I, 1 ), LDY, A( I, I ), LDA, ZERO, X( 1, I ), 1 );
               sgemv('No transpose', M-I, I-1, -ONE, A( I+1, 1 ), LDA, X( 1, I ), 1, ONE, X( I+1, I ), 1 );
               sgemv('No transpose', I-1, N-I+1, ONE, A( 1, I ), LDA, A( I, I ), LDA, ZERO, X( 1, I ), 1 );
               sgemv('No transpose', M-I, I-1, -ONE, X( I+1, 1 ), LDX, X( 1, I ), 1, ONE, X( I+1, I ), 1 );
               sscal(M-I, TAUP( I ), X( I+1, I ), 1 );

               // Update A(i+1:m,i)

               sgemv('No transpose', M-I, I-1, -ONE, A( I+1, 1 ), LDA, Y( I, 1 ), LDY, ONE, A( I+1, I ), 1 );
               sgemv('No transpose', M-I, I, -ONE, X( I+1, 1 ), LDX, A( 1, I ), 1, ONE, A( I+1, I ), 1 );

               // Generate reflection Q(i) to annihilate A(i+2:m,i)

               slarfg(M-I, A( I+1, I ), A( min( I+2, M ), I ), 1, TAUQ( I ) );
               E[I] = A( I+1, I );
               A[I+1][I] = ONE;

               // Compute Y(i+1:n,i)

               sgemv('Transpose', M-I, N-I, ONE, A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, Y( I+1, I ), 1 );
               sgemv('Transpose', M-I, I-1, ONE, A( I+1, 1 ), LDA, A( I+1, I ), 1, ZERO, Y( 1, I ), 1 );
               sgemv('No transpose', N-I, I-1, -ONE, Y( I+1, 1 ), LDY, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 );
               sgemv('Transpose', M-I, I, ONE, X( I+1, 1 ), LDX, A( I+1, I ), 1, ZERO, Y( 1, I ), 1 );
               sgemv('Transpose', I, N-I, -ONE, A( 1, I+1 ), LDA, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 );
               sscal(N-I, TAUQ( I ), Y( I+1, I ), 1 );
            }
         } // 20
      }
      }
