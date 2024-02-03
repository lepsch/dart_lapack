      SUBROUTINE DGTT02( TRANS, N, NRHS, DL, D, DU, X, LDX, B, LDB, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             TRANS;
      int                LDB, LDX, N, NRHS
      DOUBLE PRECISION   RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * ), D( * ), DL( * ), DU( * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      int                J
      DOUBLE PRECISION   ANORM, BNORM, EPS, XNORM
*     ..
*     .. External Functions ..
      bool               LSAME;
      DOUBLE PRECISION   DASUM, DLAMCH, DLANGT
      EXTERNAL           LSAME, DASUM, DLAMCH, DLANGT
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAGTM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0 or NRHS = 0
*
      RESID = ZERO
      IF( N.LE.0 .OR. NRHS.EQ.0 ) RETURN
*
*     Compute the maximum over the number of right hand sides of
*        norm(B - op(A)*X) / ( norm(op(A)) * norm(X) * EPS ).
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         ANORM = DLANGT( '1', N, DL, D, DU )
      ELSE
         ANORM = DLANGT( 'I', N, DL, D, DU )
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
*     Compute B - op(A)*X and store in B.
*
      CALL DLAGTM( TRANS, N, NRHS, -ONE, DL, D, DU, X, LDX, ONE, B, LDB )
*
      DO 10 J = 1, NRHS
         BNORM = DASUM( N, B( 1, J ), 1 )
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
*     End of DGTT02
*
      END
