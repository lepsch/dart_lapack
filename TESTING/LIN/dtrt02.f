      SUBROUTINE DTRT02( UPLO, TRANS, DIAG, N, NRHS, A, LDA, X, LDX, B, LDB, WORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            LDA, LDB, LDX, N, NRHS
      DOUBLE PRECISION   RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J
      DOUBLE PRECISION   ANORM, BNORM, EPS, XNORM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DASUM, DLAMCH, DLANTR
      EXTERNAL           LSAME, DASUM, DLAMCH, DLANTR
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DTRMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0 or NRHS = 0
*
      IF( N.LE.0 .OR. NRHS.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
*     Compute the 1-norm of op(A).
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         ANORM = DLANTR( '1', UPLO, DIAG, N, N, A, LDA, WORK )
      ELSE
         ANORM = DLANTR( 'I', UPLO, DIAG, N, N, A, LDA, WORK )
      END IF
*
*     Exit with RESID = 1/EPS if ANORM = 0.
*
      EPS = DLAMCH( 'Epsilon' )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
*
*     Compute the maximum over the number of right hand sides of
*        norm(op(A)*X - B) / ( norm(op(A)) * norm(X) * EPS )
*
      RESID = ZERO
      DO 10 J = 1, NRHS
         CALL DCOPY( N, X( 1, J ), 1, WORK, 1 )
         CALL DTRMV( UPLO, TRANS, DIAG, N, A, LDA, WORK, 1 )
         CALL DAXPY( N, -ONE, B( 1, J ), 1, WORK, 1 )
         BNORM = DASUM( N, WORK, 1 )
         XNORM = DASUM( N, X( 1, J ), 1 )
         IF( XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         END IF
   10 CONTINUE
*
      RETURN
*
*     End of DTRT02
*
      END
