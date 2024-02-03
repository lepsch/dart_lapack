      SUBROUTINE DLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                K, LDA, LDT, LDY, N, NB;
      // ..
      // .. Array Arguments ..
      double            A( LDA, * ), T( LDT, NB ), TAU( NB ), Y( LDY, NB );
      // ..

// =====================================================================

      // .. Parameters ..
      double            ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      double            EI;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DGEMM, DGEMV, DLACPY, DLARFG, DSCAL, DTRMM, DTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if (N <= 1) return;

      for (I = 1; I <= NB; I++) { // 10
         if ( I > 1 ) {

            // Update A(K+1:N,I)

            // Update I-th column of A - Y * V**T

            dgemv('NO TRANSPOSE', N-K, I-1, -ONE, Y(K+1,1), LDY, A( K+I-1, 1 ), LDA, ONE, A( K+1, I ), 1 );

            // Apply I - V * T**T * V**T to this column (call it b) from the
            // left, using the last column of T as workspace

            // Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
                     // ( V2 )             ( b2 )

            // where V1 is unit lower triangular

            // w := V1**T * b1

            dcopy(I-1, A( K+1, I ), 1, T( 1, NB ), 1 );
            dtrmv('Lower', 'Transpose', 'UNIT', I-1, A( K+1, 1 ), LDA, T( 1, NB ), 1 );

            // w := w + V2**T * b2

            dgemv('Transpose', N-K-I+1, I-1, ONE, A( K+I, 1 ), LDA, A( K+I, I ), 1, ONE, T( 1, NB ), 1 );

            // w := T**T * w

            dtrmv('Upper', 'Transpose', 'NON-UNIT', I-1, T, LDT, T( 1, NB ), 1 );

            // b2 := b2 - V2*w

            dgemv('NO TRANSPOSE', N-K-I+1, I-1, -ONE, A( K+I, 1 ), LDA, T( 1, NB ), 1, ONE, A( K+I, I ), 1 );

            // b1 := b1 - V1*w

            dtrmv('Lower', 'NO TRANSPOSE', 'UNIT', I-1, A( K+1, 1 ), LDA, T( 1, NB ), 1 );
            daxpy(I-1, -ONE, T( 1, NB ), 1, A( K+1, I ), 1 );

            A( K+I-1, I-1 ) = EI;
         }

         // Generate the elementary reflector H(I) to annihilate
         // A(K+I+1:N,I)

         dlarfg(N-K-I+1, A( K+I, I ), A( MIN( K+I+1, N ), I ), 1, TAU( I ) );
         EI = A( K+I, I );
         A( K+I, I ) = ONE;

         // Compute  Y(K+1:N,I)

         dgemv('NO TRANSPOSE', N-K, N-K-I+1, ONE, A( K+1, I+1 ), LDA, A( K+I, I ), 1, ZERO, Y( K+1, I ), 1 );
         dgemv('Transpose', N-K-I+1, I-1, ONE, A( K+I, 1 ), LDA, A( K+I, I ), 1, ZERO, T( 1, I ), 1 );
         dgemv('NO TRANSPOSE', N-K, I-1, -ONE, Y( K+1, 1 ), LDY, T( 1, I ), 1, ONE, Y( K+1, I ), 1 );
         dscal(N-K, TAU( I ), Y( K+1, I ), 1 );

         // Compute T(1:I,I)

         dscal(I-1, -TAU( I ), T( 1, I ), 1 );
         dtrmv('Upper', 'No Transpose', 'NON-UNIT', I-1, T, LDT, T( 1, I ), 1 );
         T( I, I ) = TAU( I );

      } // 10
      A( K+NB, NB ) = EI;

      // Compute Y(1:K,1:NB)

      dlacpy('ALL', K, NB, A( 1, 2 ), LDA, Y, LDY );
      dtrmm('RIGHT', 'Lower', 'NO TRANSPOSE', 'UNIT', K, NB, ONE, A( K+1, 1 ), LDA, Y, LDY )       IF( N > K+NB ) CALL DGEMM( 'NO TRANSPOSE', 'NO TRANSPOSE', K, NB, N-K-NB, ONE, A( 1, 2+NB ), LDA, A( K+1+NB, 1 ), LDA, ONE, Y, LDY );
      dtrmm('RIGHT', 'Upper', 'NO TRANSPOSE', 'NON-UNIT', K, NB, ONE, T, LDT, Y, LDY );

      return;

      // End of DLAHR2

      }
