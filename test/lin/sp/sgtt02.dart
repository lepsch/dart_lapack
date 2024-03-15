      void sgtt02(final int TRANS, final int N, final int NRHS, final int DL, final int D, final int DU, final Matrix<double> X_, final int LDX, final Matrix<double> B_, final int LDB, final int RESID,) {
  final X = X_.dim();
  final B = B_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS;
      int                LDB, LDX, N, NRHS;
      double               RESID;
      double               B( LDB, * ), D( * ), DL( * ), DU( * ), X( LDX, * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                J;
      double               ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SASUM, SLAMCH, SLANGT;
      // EXTERNAL lsame, SASUM, SLAMCH, SLANGT
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLAGTM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Quick exit if N = 0 or NRHS = 0

      RESID = ZERO;
      if (N <= 0 || NRHS == 0) return;

      // Compute the maximum over the number of right hand sides of
      //    norm(B - op(A)*X) / ( norm(op(A)) * norm(X) * EPS ).

      if ( lsame( TRANS, 'N' ) ) {
         ANORM = SLANGT( '1', N, DL, D, DU );
      } else {
         ANORM = SLANGT( 'I', N, DL, D, DU );
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Compute B - op(A)*X and store in B.

      slagtm(TRANS, N, NRHS, -ONE, DL, D, DU, X, LDX, ONE, B, LDB );

      for (J = 1; J <= NRHS; J++) { // 10
         BNORM = SASUM( N, B( 1, J ), 1 );
         XNORM = SASUM( N, X( 1, J ), 1 );
         if ( XNORM <= ZERO ) {
            RESID = ONE / EPS;
         } else {
            RESID = max( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS );
         }
      } // 10

      }