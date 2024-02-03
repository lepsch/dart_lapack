      SUBROUTINE CLATRZ( M, N, L, A, LDA, TAU, WORK )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                L, LDA, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ZERO
      const              ZERO = ( 0.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I;
      COMPLEX            ALPHA
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACGV, CLARFG, CLARZ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      IF( M.EQ.0 ) THEN
         RETURN
      ELSE IF( M.EQ.N ) THEN
         DO 10 I = 1, N
            TAU( I ) = ZERO
   10    CONTINUE
         RETURN
      END IF

      DO 20 I = M, 1, -1

         // Generate elementary reflector H(i) to annihilate
         // [ A(i,i) A(i,n-l+1:n) ]

         CALL CLACGV( L, A( I, N-L+1 ), LDA )
         ALPHA = CONJG( A( I, I ) )
         CALL CLARFG( L+1, ALPHA, A( I, N-L+1 ), LDA, TAU( I ) )
         TAU( I ) = CONJG( TAU( I ) )

         // Apply H(i) to A(1:i-1,i:n) from the right

         CALL CLARZ( 'Right', I-1, N-I+1, L, A( I, N-L+1 ), LDA, CONJG( TAU( I ) ), A( 1, I ), LDA, WORK )
         A( I, I ) = CONJG( ALPHA )

   20 CONTINUE

      RETURN

      // End of CLATRZ

      }
