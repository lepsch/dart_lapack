      SUBROUTINE ZGET03( N, A, LDA, AINV, LDAINV, WORK, LDWORK, RWORK, RCOND, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                LDA, LDAINV, LDWORK, N;
      double             RCOND, RESID;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), AINV( LDAINV, * ), WORK( LDWORK, * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) )
      // ..
      // .. Local Scalars ..
      int                I;
      double             AINVNM, ANORM, EPS;
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANGE;
      // EXTERNAL DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..
*
      // Quick exit if N = 0.
*
      IF( N.LE.0 ) THEN
         RCOND = ONE
         RESID = ZERO
         RETURN
      END IF
*
      // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.
*
      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANGE( '1', N, N, A, LDA, RWORK )
      AINVNM = ZLANGE( '1', N, N, AINV, LDAINV, RWORK )
      IF( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
         RCOND = ZERO
         RESID = ONE / EPS
         RETURN
      END IF
      RCOND = ( ONE / ANORM ) / AINVNM
*
      // Compute I - A * AINV
*
      CALL ZGEMM( 'No transpose', 'No transpose', N, N, N, -CONE, AINV, LDAINV, A, LDA, CZERO, WORK, LDWORK )
      DO 10 I = 1, N
         WORK( I, I ) = CONE + WORK( I, I )
   10 CONTINUE
*
      // Compute norm(I - AINV*A) / (N * norm(A) * norm(AINV) * EPS)
*
      RESID = ZLANGE( '1', N, N, WORK, LDWORK, RWORK )
*
      RESID = ( ( RESID*RCOND ) / EPS ) / DBLE( N )
*
      RETURN
*
      // End of ZGET03
*
      END
