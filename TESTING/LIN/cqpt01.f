      REAL             FUNCTION CQPT01( M, N, K, A, AF, LDA, TAU, JPVT, WORK, LWORK )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      int                JPVT( * )
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
      int                I, INFO, J
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
      EXTERNAL           CAXPY, CCOPY, CUNMQR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
      CQPT01 = ZERO
*
*     Test if there is enough workspace
*
      IF( LWORK.LT.M*N+N ) THEN
         CALL XERBLA( 'CQPT01', 10 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 ) RETURN
*
      NORMA = CLANGE( 'One-norm', M, N, A, LDA, RWORK )
*
      DO J = 1, K
         DO I = 1, MIN( J, M )
            WORK( ( J-1 )*M+I ) = AF( I, J )
         END DO
         DO I = J + 1, M
            WORK( ( J-1 )*M+I ) = ZERO
         END DO
      END DO
      DO J = K + 1, N
         CALL CCOPY( M, AF( 1, J ), 1, WORK( ( J-1 )*M+1 ), 1 )
      END DO
*
      CALL CUNMQR( 'Left', 'No transpose', M, N, K, AF, LDA, TAU, WORK, M, WORK( M*N+1 ), LWORK-M*N, INFO )
*
      DO J = 1, N
*
*        Compare i-th column of QR and jpvt(i)-th column of A
*
         CALL CAXPY( M, CMPLX( -ONE ), A( 1, JPVT( J ) ), 1, WORK( ( J-1 )*M+1 ), 1 )
      END DO
*
      CQPT01 = CLANGE( 'One-norm', M, N, WORK, M, RWORK ) / ( REAL( MAX( M, N ) )*SLAMCH( 'Epsilon' ) )       IF( NORMA.NE.ZERO ) CQPT01 = CQPT01 / NORMA
*
      RETURN
*
*     End of CQPT01
*
      END
