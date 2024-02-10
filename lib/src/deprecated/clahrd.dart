      void clahrd(final int N, final int K, final int NB, final Matrix<double> A, final int LDA, final int TAU, final Matrix<double> T, final int LDT, final int Y, final int LDY) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                K, LDA, LDT, LDY, N, NB;
      Complex            A( LDA, * ), T( LDT, NB ), TAU( NB ), Y( LDY, NB );
      // ..

      Complex            ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      int                I;
      Complex            EI;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CGEMV, CLACGV, CLARFG, CSCAL, CTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN

      // Quick return if possible

      if (N <= 1) return;

      for (I = 1; I <= NB; I++) { // 10
         if ( I > 1 ) {

            // Update A(1:n,i)

            // Compute i-th column of A - Y * V**H

            clacgv(I-1, A( K+I-1, 1 ), LDA );
            cgemv('No transpose', N, I-1, -ONE, Y, LDY, A( K+I-1, 1 ), LDA, ONE, A( 1, I ), 1 );
            clacgv(I-1, A( K+I-1, 1 ), LDA );

            // Apply I - V * T**H * V**H to this column (call it b) from the
            // left, using the last column of T as workspace

            // Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
            //          ( V2 )             ( b2 )

            // where V1 is unit lower triangular

            // w := V1**H * b1

            ccopy(I-1, A( K+1, I ), 1, T( 1, NB ), 1 );
            ctrmv('Lower', 'Conjugate transpose', 'Unit', I-1, A( K+1, 1 ), LDA, T( 1, NB ), 1 );

            // w := w + V2**H *b2

            cgemv('Conjugate transpose', N-K-I+1, I-1, ONE, A( K+I, 1 ), LDA, A( K+I, I ), 1, ONE, T( 1, NB ), 1 );

            // w := T**H *w

            ctrmv('Upper', 'Conjugate transpose', 'Non-unit', I-1, T, LDT, T( 1, NB ), 1 );

            // b2 := b2 - V2*w

            cgemv('No transpose', N-K-I+1, I-1, -ONE, A( K+I, 1 ), LDA, T( 1, NB ), 1, ONE, A( K+I, I ), 1 );

            // b1 := b1 - V1*w

            ctrmv('Lower', 'No transpose', 'Unit', I-1, A( K+1, 1 ), LDA, T( 1, NB ), 1 );
            caxpy(I-1, -ONE, T( 1, NB ), 1, A( K+1, I ), 1 );

            A[K+I-1][I-1] = EI;
         }

         // Generate the elementary reflector H(i) to annihilate
         // A(k+i+1:n,i)

         EI = A( K+I, I );
         clarfg(N-K-I+1, EI, A( min( K+I+1, N ), I ), 1, TAU( I ) );
         A[K+I][I] = ONE;

         // Compute  Y(1:n,i)

         cgemv('No transpose', N, N-K-I+1, ONE, A( 1, I+1 ), LDA, A( K+I, I ), 1, ZERO, Y( 1, I ), 1 );
         cgemv('Conjugate transpose', N-K-I+1, I-1, ONE, A( K+I, 1 ), LDA, A( K+I, I ), 1, ZERO, T( 1, I ), 1 );
         cgemv('No transpose', N, I-1, -ONE, Y, LDY, T( 1, I ), 1, ONE, Y( 1, I ), 1 );
         cscal(N, TAU( I ), Y( 1, I ), 1 );

         // Compute T(1:i,i)

         cscal(I-1, -TAU( I ), T( 1, I ), 1 );
         ctrmv('Upper', 'No transpose', 'Non-unit', I-1, T, LDT, T( 1, I ), 1 );
         T[I][I] = TAU( I );

      } // 10
      A[K+NB][NB] = EI;

      }
