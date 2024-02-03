      SUBROUTINE ZTRT02( UPLO, TRANS, DIAG, N, NRHS, A, LDA, X, LDX, B, LDB, WORK, RWORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                LDA, LDB, LDX, N, NRHS
      double             RESID;
*     ..
*     .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      int                J
      double             ANORM, BNORM, EPS, XNORM;
*     ..
*     .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DZASUM, ZLANTR;
      EXTERNAL           LSAME, DLAMCH, DZASUM, ZLANTR
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZAXPY, ZCOPY, ZTRMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, MAX
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
         ANORM = ZLANTR( '1', UPLO, DIAG, N, N, A, LDA, RWORK )
      ELSE
         ANORM = ZLANTR( 'I', UPLO, DIAG, N, N, A, LDA, RWORK )
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
         CALL ZCOPY( N, X( 1, J ), 1, WORK, 1 )
         CALL ZTRMV( UPLO, TRANS, DIAG, N, A, LDA, WORK, 1 )
         CALL ZAXPY( N, DCMPLX( -ONE ), B( 1, J ), 1, WORK, 1 )
         BNORM = DZASUM( N, WORK, 1 )
         XNORM = DZASUM( N, X( 1, J ), 1 )
         IF( XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         END IF
   10 CONTINUE
*
      RETURN
*
*     End of ZTRT02
*
      END
