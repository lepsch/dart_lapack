      SUBROUTINE SGTT02( TRANS, N, NRHS, DL, D, DU, X, LDX, B, LDB, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                LDB, LDX, N, NRHS;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      REAL               B( LDB, * ), D( * ), DL( * ), DU( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                J;
      REAL               ANORM, BNORM, EPS, XNORM
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SASUM, SLAMCH, SLANGT
      // EXTERNAL LSAME, SASUM, SLAMCH, SLANGT
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLAGTM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0

      RESID = ZERO
      if (N <= 0 || NRHS == 0) RETURN;

      // Compute the maximum over the number of right hand sides of
         // norm(B - op(A)*X) / ( norm(op(A)) * norm(X) * EPS ).

      if ( LSAME( TRANS, 'N' ) ) {
         ANORM = SLANGT( '1', N, DL, D, DU )
      } else {
         ANORM = SLANGT( 'I', N, DL, D, DU )
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' )
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS
         RETURN
      }

      // Compute B - op(A)*X and store in B.

      slagtm(TRANS, N, NRHS, -ONE, DL, D, DU, X, LDX, ONE, B, LDB );

      for (J = 1; J <= NRHS; J++) { // 10
         BNORM = SASUM( N, B( 1, J ), 1 )
         XNORM = SASUM( N, X( 1, J ), 1 )
         if ( XNORM <= ZERO ) {
            RESID = ONE / EPS
         } else {
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         }
      } // 10

      RETURN

      // End of SGTT02

      }
