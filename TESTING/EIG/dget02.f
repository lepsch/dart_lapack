      SUBROUTINE DGET02( TRANS, M, N, NRHS, A, LDA, X, LDX, B, LDB, RWORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            LDA, LDB, LDX, M, N, NRHS
      DOUBLE PRECISION   RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), RWORK( * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, N1, N2
      DOUBLE PRECISION   ANORM, BNORM, EPS, XNORM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DASUM, DLAMCH, DLANGE
      EXTERNAL           LSAME, DASUM, DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Quick exit if M = 0 or N = 0 or NRHS = 0
*
      IF( M.LE.0 .OR. N.LE.0 .OR. NRHS.EQ.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
      IF( LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' ) ) THEN
         N1 = N
         N2 = M
      ELSE
         N1 = M
         N2 = N
      END IF
*
*     Exit with RESID = 1/EPS if ANORM = 0.
*
      EPS = DLAMCH( 'Epsilon' )
      IF( LSAME( TRANS, 'N' ) ) THEN
         ANORM = DLANGE( '1', M, N, A, LDA, RWORK )
      ELSE
         ANORM = DLANGE( 'I', M, N, A, LDA, RWORK )
      END IF
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
*
*     Compute B - op(A)*X and store in B.
*
      CALL DGEMM( TRANS, 'No transpose', N1, NRHS, N2, -ONE, A, LDA, X, LDX, ONE, B, LDB )
*
*     Compute the maximum over the number of right hand sides of
*        norm(B - op(A)*X) / ( norm(op(A)) * norm(X) * EPS ) .
*
      RESID = ZERO
      DO 10 J = 1, NRHS
         BNORM = DASUM( N1, B( 1, J ), 1 )
         XNORM = DASUM( N2, X( 1, J ), 1 )
         IF( XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         END IF
   10 CONTINUE
*
      RETURN
*
*     End of DGET02
*
      END
