      REAL             FUNCTION CRZT01( M, N, A, AF, LDA, TAU, WORK, LWORK )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), AF( LDA, * ), TAU( * ), WORK( LWORK )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J
      REAL               NORMA
*     ..
*     .. Local Arrays ..
      REAL               RWORK( 1 )
*     ..
*     .. External Functions ..
      REAL               CLANGE, SLAMCH
      EXTERNAL           CLANGE, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CAXPY, CLASET, CUNMRZ, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, MAX, REAL
*     ..
*     .. Executable Statements ..
*
      CRZT01 = ZERO
*
      IF( LWORK.LT.M*N+M ) THEN
         CALL XERBLA( 'CRZT01', 8 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 ) RETURN
*
      NORMA = CLANGE( 'One-norm', M, N, A, LDA, RWORK )
*
*     Copy upper triangle R
*
      CALL CLASET( 'Full', M, N, CMPLX( ZERO ), CMPLX( ZERO ), WORK, M )
      DO 20 J = 1, M
         DO 10 I = 1, J
            WORK( ( J-1 )*M+I ) = AF( I, J )
   10    CONTINUE
   20 CONTINUE
*
*     R = R * P(1) * ... *P(m)
*
      CALL CUNMRZ( 'Right', 'No transpose', M, N, M, N-M, AF, LDA, TAU, WORK, M, WORK( M*N+1 ), LWORK-M*N, INFO )
*
*     R = R - A
*
      DO 30 I = 1, N
         CALL CAXPY( M, CMPLX( -ONE ), A( 1, I ), 1, WORK( ( I-1 )*M+1 ), 1 )
   30 CONTINUE
*
      CRZT01 = CLANGE( 'One-norm', M, N, WORK, M, RWORK )
*
      CRZT01 = CRZT01 / ( SLAMCH( 'Epsilon' )*REAL( MAX( M, N ) ) )
      IF( NORMA.NE.ZERO ) CRZT01 = CRZT01 / NORMA
*
      RETURN
*
*     End of CRZT01
*
      END
