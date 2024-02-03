      SUBROUTINE ZGTT02( TRANS, N, NRHS, DL, D, DU, X, LDX, B, LDB, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                LDB, LDX, N, NRHS;
      double             RESID;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         B( LDB, * ), D( * ), DL( * ), DU( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      // ..
      // .. Local Scalars ..
      int                J;
      double             ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DZASUM, ZLANGT;
      // EXTERNAL LSAME, DLAMCH, DZASUM, ZLANGT
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLAGTM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0

      RESID = ZERO
      IF( N.LE.0 .OR. NRHS.EQ.0 ) RETURN

      // Compute the maximum over the number of right hand sides of
         // norm(B - op(A)*X) / ( norm(op(A)) * norm(X) * EPS ).

      IF( LSAME( TRANS, 'N' ) ) THEN
         ANORM = ZLANGT( '1', N, DL, D, DU )
      ELSE
         ANORM = ZLANGT( 'I', N, DL, D, DU )
      END IF

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = DLAMCH( 'Epsilon' )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF

      // Compute B - op(A)*X and store in B.

      CALL ZLAGTM( TRANS, N, NRHS, -ONE, DL, D, DU, X, LDX, ONE, B, LDB )

      DO 10 J = 1, NRHS
         BNORM = DZASUM( N, B( 1, J ), 1 )
         XNORM = DZASUM( N, X( 1, J ), 1 )
         IF( XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         END IF
   10 CONTINUE

      RETURN

      // End of ZGTT02

      }
