      SUBROUTINE ZLATRZ( M, N, L, A, LDA, TAU, WORK )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                L, LDA, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ZERO
      const              ZERO = ( 0.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I;
      COMPLEX*16         ALPHA
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLACGV, ZLARFG, ZLARZ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( M.EQ.0 ) {
         RETURN
      } else if ( M.EQ.N ) {
         DO 10 I = 1, N
            TAU( I ) = ZERO
   10    CONTINUE
         RETURN
      }

      DO 20 I = M, 1, -1

         // Generate elementary reflector H(i) to annihilate
         // [ A(i,i) A(i,n-l+1:n) ]

         CALL ZLACGV( L, A( I, N-L+1 ), LDA )
         ALPHA = DCONJG( A( I, I ) )
         CALL ZLARFG( L+1, ALPHA, A( I, N-L+1 ), LDA, TAU( I ) )
         TAU( I ) = DCONJG( TAU( I ) )

         // Apply H(i) to A(1:i-1,i:n) from the right

         CALL ZLARZ( 'Right', I-1, N-I+1, L, A( I, N-L+1 ), LDA, DCONJG( TAU( I ) ), A( 1, I ), LDA, WORK )
         A( I, I ) = DCONJG( ALPHA )

   20 CONTINUE

      RETURN

      // End of ZLATRZ

      }
