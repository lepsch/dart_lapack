      SUBROUTINE SGTT02( TRANS, N, NRHS, DL, D, DU, X, LDX, B, LDB, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             TRANS;
      int                LDB, LDX, N, NRHS
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               B( LDB, * ), D( * ), DL( * ), DU( * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      int                J
      REAL               ANORM, BNORM, EPS, XNORM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SASUM, SLAMCH, SLANGT
      EXTERNAL           LSAME, SASUM, SLAMCH, SLANGT
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLAGTM
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
         ANORM = SLANGT( '1', N, DL, D, DU )
      ELSE
         ANORM = SLANGT( 'I', N, DL, D, DU )
      END IF
*
*     Exit with RESID = 1/EPS if ANORM = 0.
*
      EPS = SLAMCH( 'Epsilon' )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
*
*     Compute B - op(A)*X and store in B.
*
      CALL SLAGTM( TRANS, N, NRHS, -ONE, DL, D, DU, X, LDX, ONE, B, LDB )
*
      DO 10 J = 1, NRHS
         BNORM = SASUM( N, B( 1, J ), 1 )
         XNORM = SASUM( N, X( 1, J ), 1 )
         IF( XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         END IF
   10 CONTINUE
*
      RETURN
*
*     End of SGTT02
*
      END
