      void slahrd(N, K, NB, A, LDA, TAU, T, LDT, Y, LDY ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                K, LDA, LDT, LDY, N, NB;
      // ..
      // .. Array Arguments ..
      double               A( LDA, * ), T( LDT, NB ), TAU( NB ), Y( LDY, NB );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      double               EI;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, SGEMV, SLARFG, SSCAL, STRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if (N <= 1) return;

      for (I = 1; I <= NB; I++) { // 10
         if ( I > 1 ) {

            // Update A(1:n,i)

            // Compute i-th column of A - Y * V**T

            sgemv('No transpose', N, I-1, -ONE, Y, LDY, A( K+I-1, 1 ), LDA, ONE, A( 1, I ), 1 );

            // Apply I - V * T**T * V**T to this column (call it b) from the
            // left, using the last column of T as workspace

            // Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
                     // ( V2 )             ( b2 )

            // where V1 is unit lower triangular

            // w := V1**T * b1

            scopy(I-1, A( K+1, I ), 1, T( 1, NB ), 1 );
            strmv('Lower', 'Transpose', 'Unit', I-1, A( K+1, 1 ), LDA, T( 1, NB ), 1 );

            // w := w + V2**T *b2

            sgemv('Transpose', N-K-I+1, I-1, ONE, A( K+I, 1 ), LDA, A( K+I, I ), 1, ONE, T( 1, NB ), 1 );

            // w := T**T *w

            strmv('Upper', 'Transpose', 'Non-unit', I-1, T, LDT, T( 1, NB ), 1 );

            // b2 := b2 - V2*w

            sgemv('No transpose', N-K-I+1, I-1, -ONE, A( K+I, 1 ), LDA, T( 1, NB ), 1, ONE, A( K+I, I ), 1 );

            // b1 := b1 - V1*w

            strmv('Lower', 'No transpose', 'Unit', I-1, A( K+1, 1 ), LDA, T( 1, NB ), 1 );
            saxpy(I-1, -ONE, T( 1, NB ), 1, A( K+1, I ), 1 );

            A[K+I-1, I-1] = EI;
         }

         // Generate the elementary reflector H(i) to annihilate
         // A(k+i+1:n,i)

         slarfg(N-K-I+1, A( K+I, I ), A( min( K+I+1, N ), I ), 1, TAU( I ) );
         EI = A( K+I, I );
         A[K+I, I] = ONE;

         // Compute  Y(1:n,i)

         sgemv('No transpose', N, N-K-I+1, ONE, A( 1, I+1 ), LDA, A( K+I, I ), 1, ZERO, Y( 1, I ), 1 );
         sgemv('Transpose', N-K-I+1, I-1, ONE, A( K+I, 1 ), LDA, A( K+I, I ), 1, ZERO, T( 1, I ), 1 );
         sgemv('No transpose', N, I-1, -ONE, Y, LDY, T( 1, I ), 1, ONE, Y( 1, I ), 1 );
         sscal(N, TAU( I ), Y( 1, I ), 1 );

         // Compute T(1:i,i)

         sscal(I-1, -TAU( I ), T( 1, I ), 1 );
         strmv('Upper', 'No transpose', 'Non-unit', I-1, T, LDT, T( 1, I ), 1 );
         T[I, I] = TAU( I );

      } // 10
      A[K+NB, NB] = EI;

      return;
      }
