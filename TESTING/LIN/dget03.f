      SUBROUTINE DGET03( N, A, LDA, AINV, LDAINV, WORK, LDWORK, RWORK, RCOND, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDA, LDAINV, LDWORK, N;
      double             RCOND, RESID;
*     ..
*     .. Array Arguments ..
      double             A( LDA, * ), AINV( LDAINV, * ), RWORK( * ), WORK( LDWORK, * );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      int                I;
      double             AINVNM, ANORM, EPS;
*     ..
*     .. External Functions ..
      double             DLAMCH, DLANGE;
      EXTERNAL           DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC DBLE
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
      EPS = DLAMCH( 'Epsilon' )
      ANORM = DLANGE( '1', N, N, A, LDA, RWORK )
      AINVNM = DLANGE( '1', N, N, AINV, LDAINV, RWORK )
      IF( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
         RCOND = ZERO
         RESID = ONE / EPS
         RETURN
      END IF
      RCOND = ( ONE / ANORM ) / AINVNM
*
*     Compute I - A * AINV
*
      CALL DGEMM( 'No transpose', 'No transpose', N, N, N, -ONE, AINV, LDAINV, A, LDA, ZERO, WORK, LDWORK )
      DO 10 I = 1, N
         WORK( I, I ) = ONE + WORK( I, I )
   10 CONTINUE
*
*     Compute norm(I - AINV*A) / (N * norm(A) * norm(AINV) * EPS)
*
      RESID = DLANGE( '1', N, N, WORK, LDWORK, RWORK )
*
      RESID = ( ( RESID*RCOND ) / EPS ) / DBLE( N )
*
      RETURN
*
*     End of DGET03
*
      END
