      SUBROUTINE CGET03( N, A, LDA, AINV, LDAINV, WORK, LDWORK, RWORK, RCOND, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDAINV, LDWORK, N
      REAL               RCOND, RESID
*     ..
*     .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), AINV( LDAINV, * ), WORK( LDWORK, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      REAL               AINVNM, ANORM, EPS
*     ..
*     .. External Functions ..
      REAL               CLANGE, SLAMCH
      EXTERNAL           CLANGE, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0.
*
      IF( N.LE.0 ) THEN
         RCOND = ONE
         RESID = ZERO
         RETURN
      END IF
*
*     Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.
*
      EPS = SLAMCH( 'Epsilon' )
      ANORM = CLANGE( '1', N, N, A, LDA, RWORK )
      AINVNM = CLANGE( '1', N, N, AINV, LDAINV, RWORK )
      IF( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
         RCOND = ZERO
         RESID = ONE / EPS
         RETURN
      END IF
      RCOND = ( ONE/ANORM ) / AINVNM
*
*     Compute I - A * AINV
*
      CALL CGEMM( 'No transpose', 'No transpose', N, N, N, -CONE, AINV, LDAINV, A, LDA, CZERO, WORK, LDWORK )
      DO 10 I = 1, N
         WORK( I, I ) = CONE + WORK( I, I )
   10 CONTINUE
*
*     Compute norm(I - AINV*A) / (N * norm(A) * norm(AINV) * EPS)
*
      RESID = CLANGE( '1', N, N, WORK, LDWORK, RWORK )
*
      RESID = ( ( RESID*RCOND )/EPS ) / REAL( N )
*
      RETURN
*
*     End of CGET03
*
      END
