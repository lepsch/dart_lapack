      SUBROUTINE DPOT06( UPLO, N, NRHS, A, LDA, X, LDX, B, LDB, RWORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDB, LDX, N, NRHS
      DOUBLE PRECISION   RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), RWORK( * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, NEGONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      PARAMETER          ( NEGONE = -1.0D+0 )
*     ..
*     .. Local Scalars ..
      int                IFAIL, J
      DOUBLE PRECISION   ANORM, BNORM, EPS, XNORM
*     ..
*     .. External Functions ..
      int                IDAMAX
      DOUBLE PRECISION   DLAMCH, DLANSY
      EXTERNAL           IDAMAX, DLAMCH, DLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSYMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, ABS
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0 or NRHS = 0
*
      IF( N.LE.0 .OR. NRHS.EQ.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
*     Exit with RESID = 1/EPS if ANORM = 0.
*
      EPS = DLAMCH( 'Epsilon' )
      ANORM = DLANSY( 'I', UPLO, N, A, LDA, RWORK )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
*
*     Compute  B - A*X  and store in B.
      IFAIL=0
*
      CALL DSYMM( 'Left', UPLO, N, NRHS, NEGONE, A, LDA, X, LDX, ONE, B, LDB )
*
*     Compute the maximum over the number of right hand sides of
*        norm(B - A*X) / ( norm(A) * norm(X) * EPS ) .
*
      RESID = ZERO
      DO 10 J = 1, NRHS
         BNORM = ABS(B(IDAMAX( N, B( 1, J ), 1 ),J))
         XNORM = ABS(X(IDAMAX( N, X( 1, J ), 1 ),J))
         IF( XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         END IF
   10 CONTINUE
*
      RETURN
*
*     End of DPOT06
*
      END
