      void zqrt16(TRANS, M, N, NRHS, final Matrix<double> A, final int LDA, final Matrix<double> X, final int LDX, final Matrix<double> B, final int LDB, RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS;
      int                LDA, LDB, LDX, M, N, NRHS;
      double             RESID;
      double             RWORK( * );
      Complex         A( LDA, * ), B( LDB, * ), X( LDX, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      int                J, N1, N2;
      double             ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, DZASUM, ZLANGE;
      // EXTERNAL lsame, DLAMCH, DZASUM, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Quick exit if M = 0 or N = 0 or NRHS = 0

      if ( M <= 0 || N <= 0 || NRHS == 0 ) {
         RESID = ZERO;
         return;
      }

      if ( lsame( TRANS, 'T' ) || lsame( TRANS, 'C' ) ) {
         ANORM = ZLANGE( 'I', M, N, A, LDA, RWORK );
         N1 = N;
         N2 = M;
      } else {
         ANORM = ZLANGE( '1', M, N, A, LDA, RWORK );
         N1 = M;
         N2 = N;
      }

      EPS = dlamch( 'Epsilon' );

      // Compute  B - A*X  (or  B - A'*X ) and store in B.

      zgemm(TRANS, 'No transpose', N1, NRHS, N2, -CONE, A, LDA, X, LDX, CONE, B, LDB );

      // Compute the maximum over the number of right hand sides of
      //    norm(B - A*X) / ( max(m,n) * norm(A) * norm(X) * EPS ) .

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 10
         BNORM = DZASUM( N1, B( 1, J ), 1 );
         XNORM = DZASUM( N2, X( 1, J ), 1 );
         if ( ANORM == ZERO && BNORM == ZERO ) {
            RESID = ZERO;
         } else if ( ANORM <= ZERO || XNORM <= ZERO ) {
            RESID = ONE / EPS;
         } else {
            RESID = max( RESID, ( ( BNORM / ANORM ) / XNORM ) / ( max( M, N )*EPS ) );
         }
      } // 10

      }
