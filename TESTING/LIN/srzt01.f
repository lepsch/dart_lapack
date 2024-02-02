      REAL             FUNCTION SRZT01( M, N, A, AF, LDA, TAU, WORK,
     $                 LWORK )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), AF( LDA, * ), TAU( * ),
     $                   WORK( LWORK )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J
      REAL               NORMA
*     ..
*     .. Local Arrays ..
      REAL               RWORK( 1 )
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLANGE
      EXTERNAL           SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SAXPY, SLASET, SORMRZ, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, REAL
*     ..
*     .. Executable Statements ..
*
      SRZT01 = ZERO
*
      IF( LWORK.LT.M*N+M ) THEN
         CALL XERBLA( 'SRZT01', 8 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      NORMA = SLANGE( 'One-norm', M, N, A, LDA, RWORK )
*
*     Copy upper triangle R
*
      CALL SLASET( 'Full', M, N, ZERO, ZERO, WORK, M )
      DO 20 J = 1, M
         DO 10 I = 1, J
            WORK( ( J-1 )*M+I ) = AF( I, J )
   10    CONTINUE
   20 CONTINUE
*
*     R = R * P(1) * ... *P(m)
*
      CALL SORMRZ( 'Right', 'No transpose', M, N, M, N-M, AF, LDA, TAU,
     $             WORK, M, WORK( M*N+1 ), LWORK-M*N, INFO )
*
*     R = R - A
*
      DO 30 I = 1, N
         CALL SAXPY( M, -ONE, A( 1, I ), 1, WORK( ( I-1 )*M+1 ), 1 )
   30 CONTINUE
*
      SRZT01 = SLANGE( 'One-norm', M, N, WORK, M, RWORK )
*
      SRZT01 = SRZT01 / ( SLAMCH( 'Epsilon' )*REAL( MAX( M, N ) ) )
      IF( NORMA.NE.ZERO )
     $   SRZT01 = SRZT01 / NORMA
*
      RETURN
*
*     End of SRZT01
*
      END