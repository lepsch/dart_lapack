      SUBROUTINE CLAHRD( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                K, LDA, LDT, LDY, N, NB;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), T( LDT, NB ), TAU( NB ), Y( LDY, NB )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ), ONE = ( 1.0E+0, 0.0E+0 ) )
      // ..
      // .. Local Scalars ..
      int                I;
      COMPLEX            EI
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CGEMV, CLACGV, CLARFG, CSCAL, CTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN
      // ..
      // .. Executable Statements ..
*
      // Quick return if possible
*
      IF( N.LE.1 ) RETURN
*
      DO 10 I = 1, NB
         IF( I.GT.1 ) THEN
*
            // Update A(1:n,i)
*
            // Compute i-th column of A - Y * V**H
*
            CALL CLACGV( I-1, A( K+I-1, 1 ), LDA )
            CALL CGEMV( 'No transpose', N, I-1, -ONE, Y, LDY, A( K+I-1, 1 ), LDA, ONE, A( 1, I ), 1 )
            CALL CLACGV( I-1, A( K+I-1, 1 ), LDA )
*
            // Apply I - V * T**H * V**H to this column (call it b) from the
            // left, using the last column of T as workspace
*
            // Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
                     // ( V2 )             ( b2 )
*
            // where V1 is unit lower triangular
*
            // w := V1**H * b1
*
            CALL CCOPY( I-1, A( K+1, I ), 1, T( 1, NB ), 1 )
            CALL CTRMV( 'Lower', 'Conjugate transpose', 'Unit', I-1, A( K+1, 1 ), LDA, T( 1, NB ), 1 )
*
            // w := w + V2**H *b2
*
            CALL CGEMV( 'Conjugate transpose', N-K-I+1, I-1, ONE, A( K+I, 1 ), LDA, A( K+I, I ), 1, ONE, T( 1, NB ), 1 )
*
            // w := T**H *w
*
            CALL CTRMV( 'Upper', 'Conjugate transpose', 'Non-unit', I-1, T, LDT, T( 1, NB ), 1 )
*
            // b2 := b2 - V2*w
*
            CALL CGEMV( 'No transpose', N-K-I+1, I-1, -ONE, A( K+I, 1 ), LDA, T( 1, NB ), 1, ONE, A( K+I, I ), 1 )
*
            // b1 := b1 - V1*w
*
            CALL CTRMV( 'Lower', 'No transpose', 'Unit', I-1, A( K+1, 1 ), LDA, T( 1, NB ), 1 )
            CALL CAXPY( I-1, -ONE, T( 1, NB ), 1, A( K+1, I ), 1 )
*
            A( K+I-1, I-1 ) = EI
         END IF
*
         // Generate the elementary reflector H(i) to annihilate
         // A(k+i+1:n,i)
*
         EI = A( K+I, I )
         CALL CLARFG( N-K-I+1, EI, A( MIN( K+I+1, N ), I ), 1, TAU( I ) )
         A( K+I, I ) = ONE
*
         // Compute  Y(1:n,i)
*
         CALL CGEMV( 'No transpose', N, N-K-I+1, ONE, A( 1, I+1 ), LDA, A( K+I, I ), 1, ZERO, Y( 1, I ), 1 )          CALL CGEMV( 'Conjugate transpose', N-K-I+1, I-1, ONE, A( K+I, 1 ), LDA, A( K+I, I ), 1, ZERO, T( 1, I ), 1 )
         CALL CGEMV( 'No transpose', N, I-1, -ONE, Y, LDY, T( 1, I ), 1, ONE, Y( 1, I ), 1 )
         CALL CSCAL( N, TAU( I ), Y( 1, I ), 1 )
*
         // Compute T(1:i,i)
*
         CALL CSCAL( I-1, -TAU( I ), T( 1, I ), 1 )
         CALL CTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T, LDT, T( 1, I ), 1 )
         T( I, I ) = TAU( I )
*
   10 CONTINUE
      A( K+NB, NB ) = EI
*
      RETURN
*
      // End of CLAHRD
*
      END
