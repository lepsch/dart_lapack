      SUBROUTINE DGTT02( TRANS, N, NRHS, DL, D, DU, X, LDX, B, LDB, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                LDB, LDX, N, NRHS;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             B( LDB, * ), D( * ), DL( * ), DU( * ), X( LDX, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                J;
      double             ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DASUM, DLAMCH, DLANGT;
      // EXTERNAL LSAME, DASUM, DLAMCH, DLANGT
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLAGTM
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

      if ( LSAME( TRANS, 'N' ) ) {
         ANORM = DLANGT( '1', N, DL, D, DU )
      } else {
         ANORM = DLANGT( 'I', N, DL, D, DU )
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = DLAMCH( 'Epsilon' )
      if ( ANORM.LE.ZERO ) {
         RESID = ONE / EPS
         RETURN
      }

      // Compute B - op(A)*X and store in B.

      dlagtm(TRANS, N, NRHS, -ONE, DL, D, DU, X, LDX, ONE, B, LDB );

      for (J = 1; J <= NRHS; J++) { // 10
         BNORM = DASUM( N, B( 1, J ), 1 )
         XNORM = DASUM( N, X( 1, J ), 1 )
         if ( XNORM.LE.ZERO ) {
            RESID = ONE / EPS
         } else {
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         }
      } // 10

      RETURN

      // End of DGTT02

      }
