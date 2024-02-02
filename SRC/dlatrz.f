      SUBROUTINE DLATRZ( M, N, L, A, LDA, TAU, WORK )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            L, LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFG, DLARZ
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
*     Quick return if possible
*
      IF( M.EQ.0 ) THEN
         RETURN
      ELSE IF( M.EQ.N ) THEN
         DO 10 I = 1, N
            TAU( I ) = ZERO
   10    CONTINUE
         RETURN
      END IF
*
      DO 20 I = M, 1, -1
*
*        Generate elementary reflector H(i) to annihilate
*        [ A(i,i) A(i,n-l+1:n) ]
*
         CALL DLARFG( L+1, A( I, I ), A( I, N-L+1 ), LDA, TAU( I ) )
*
*        Apply H(i) to A(1:i-1,i:n) from the right
*
         CALL DLARZ( 'Right', I-1, N-I+1, L, A( I, N-L+1 ), LDA,
     $               TAU( I ), A( 1, I ), LDA, WORK )
*
   20 CONTINUE
*
      RETURN
*
*     End of DLATRZ
*
      END