      SUBROUTINE CLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                K, LDA, LDT, LDY, N, NB;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), T( LDT, NB ), TAU( NB ), Y( LDY, NB )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ZERO, ONE
      const              ZERO = ( 0.0E+0, 0.0E+0 ), ONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I;
      COMPLEX            EI
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CGEMM, CGEMV, CLACPY, CLARFG, CSCAL, CTRMM, CTRMV, CLACGV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      IF( N.LE.1 ) RETURN

      for (I = 1; I <= NB; I++) { // 10
         if ( I.GT.1 ) {

            // Update A(K+1:N,I)

            // Update I-th column of A - Y * V**H

            clacgv(I-1, A( K+I-1, 1 ), LDA );
            cgemv('NO TRANSPOSE', N-K, I-1, -ONE, Y(K+1,1), LDY, A( K+I-1, 1 ), LDA, ONE, A( K+1, I ), 1 );
            clacgv(I-1, A( K+I-1, 1 ), LDA );

            // Apply I - V * T**H * V**H to this column (call it b) from the
            // left, using the last column of T as workspace

            // Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
                     // ( V2 )             ( b2 )

            // where V1 is unit lower triangular

            // w := V1**H * b1

            ccopy(I-1, A( K+1, I ), 1, T( 1, NB ), 1 );
            ctrmv('Lower', 'Conjugate transpose', 'UNIT', I-1, A( K+1, 1 ), LDA, T( 1, NB ), 1 );

            // w := w + V2**H * b2

            cgemv('Conjugate transpose', N-K-I+1, I-1, ONE, A( K+I, 1 ), LDA, A( K+I, I ), 1, ONE, T( 1, NB ), 1 );

            // w := T**H * w

            ctrmv('Upper', 'Conjugate transpose', 'NON-UNIT', I-1, T, LDT, T( 1, NB ), 1 );

            // b2 := b2 - V2*w

            cgemv('NO TRANSPOSE', N-K-I+1, I-1, -ONE, A( K+I, 1 ), LDA, T( 1, NB ), 1, ONE, A( K+I, I ), 1 );

            // b1 := b1 - V1*w

            ctrmv('Lower', 'NO TRANSPOSE', 'UNIT', I-1, A( K+1, 1 ), LDA, T( 1, NB ), 1 );
            caxpy(I-1, -ONE, T( 1, NB ), 1, A( K+1, I ), 1 );

            A( K+I-1, I-1 ) = EI
         }

         // Generate the elementary reflector H(I) to annihilate
         // A(K+I+1:N,I)

         clarfg(N-K-I+1, A( K+I, I ), A( MIN( K+I+1, N ), I ), 1, TAU( I ) );
         EI = A( K+I, I )
         A( K+I, I ) = ONE

         // Compute  Y(K+1:N,I)

         cgemv('NO TRANSPOSE', N-K, N-K-I+1, ONE, A( K+1, I+1 ), LDA, A( K+I, I ), 1, ZERO, Y( K+1, I ), 1 )          CALL CGEMV( 'Conjugate transpose', N-K-I+1, I-1, ONE, A( K+I, 1 ), LDA, A( K+I, I ), 1, ZERO, T( 1, I ), 1 )          CALL CGEMV( 'NO TRANSPOSE', N-K, I-1, -ONE, Y( K+1, 1 ), LDY, T( 1, I ), 1, ONE, Y( K+1, I ), 1 );
         cscal(N-K, TAU( I ), Y( K+1, I ), 1 );

         // Compute T(1:I,I)

         cscal(I-1, -TAU( I ), T( 1, I ), 1 );
         ctrmv('Upper', 'No Transpose', 'NON-UNIT', I-1, T, LDT, T( 1, I ), 1 );
         T( I, I ) = TAU( I )

      } // 10
      A( K+NB, NB ) = EI

      // Compute Y(1:K,1:NB)

      clacpy('ALL', K, NB, A( 1, 2 ), LDA, Y, LDY );
      ctrmm('RIGHT', 'Lower', 'NO TRANSPOSE', 'UNIT', K, NB, ONE, A( K+1, 1 ), LDA, Y, LDY )       IF( N.GT.K+NB ) CALL CGEMM( 'NO TRANSPOSE', 'NO TRANSPOSE', K, NB, N-K-NB, ONE, A( 1, 2+NB ), LDA, A( K+1+NB, 1 ), LDA, ONE, Y, LDY );
      ctrmm('RIGHT', 'Upper', 'NO TRANSPOSE', 'NON-UNIT', K, NB, ONE, T, LDT, Y, LDY );

      RETURN

      // End of CLAHR2

      }
