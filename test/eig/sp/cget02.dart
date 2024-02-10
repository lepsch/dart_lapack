      void cget02(final int TRANS, final int M, final int N, final int NRHS, final Matrix<double> A, final int LDA, final Matrix<double> X, final int LDX, final Matrix<double> B, final int LDB, final Array<double> RWORK, final int RESID) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS;
      int                LDA, LDB, LDX, M, N, NRHS;
      double               RESID;
      double               RWORK( * );
      Complex            A( LDA, * ), B( LDB, * ), X( LDX, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      int                J, N1, N2;
      double               ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANGE, SCASUM, SLAMCH;
      // EXTERNAL lsame, CLANGE, SCASUM, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Quick exit if M = 0 or N = 0 or NRHS = 0

      if ( M <= 0 || N <= 0 || NRHS == 0 ) {
         RESID = ZERO;
         return;
      }

      if ( lsame( TRANS, 'T' ) || lsame( TRANS, 'C' ) ) {
         N1 = N;
         N2 = M;
      } else {
         N1 = M;
         N2 = N;
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' );
      if ( lsame( TRANS, 'N' ) ) {
         ANORM = CLANGE( '1', M, N, A, LDA, RWORK );
      } else {
         ANORM = CLANGE( 'I', M, N, A, LDA, RWORK );
      }
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Compute B - op(A)*X and store in B.

      cgemm(TRANS, 'No transpose', N1, NRHS, N2, -CONE, A, LDA, X, LDX, CONE, B, LDB );

      // Compute the maximum over the number of right hand sides of
      //    norm(B - op(A)*X) / ( norm(op(A)) * norm(X) * EPS ) .

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 10
         BNORM = SCASUM( N1, B( 1, J ), 1 );
         XNORM = SCASUM( N2, X( 1, J ), 1 );
         if ( XNORM <= ZERO ) {
            RESID = ONE / EPS;
         } else {
            RESID = max( RESID, ( ( BNORM/ANORM )/XNORM )/EPS );
         }
      } // 10

      }
