      void dqrt16(final int TRANS, final int M, final int N, final int NRHS, final Matrix<double> A_, final int LDA, final Matrix<double> X_, final int LDX, final Matrix<double> B_, final int LDB, final Array<double> RWORK_, final int RESID,) {
  final A = A_.dim();
  final X = X_.dim();
  final B = B_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS;
      int                LDA, LDB, LDX, M, N, NRHS;
      double             RESID;
      double             A( LDA, * ), B( LDB, * ), RWORK( * ), X( LDX, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                J, N1, N2;
      double             ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DASUM, DLAMCH, DLANGE;
      // EXTERNAL lsame, DASUM, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Quick exit if M = 0 or N = 0 or NRHS = 0

      if ( M <= 0 || N <= 0 || NRHS == 0 ) {
         RESID = ZERO;
         return;
      }

      if ( lsame( TRANS, 'T' ) || lsame( TRANS, 'C' ) ) {
         ANORM = dlange( 'I', M, N, A, LDA, RWORK );
         N1 = N;
         N2 = M;
      } else {
         ANORM = dlange( '1', M, N, A, LDA, RWORK );
         N1 = M;
         N2 = N;
      }

      EPS = dlamch( 'Epsilon' );

      // Compute  B - A*X  (or  B - A'*X ) and store in B.

      dgemm(TRANS, 'No transpose', N1, NRHS, N2, -ONE, A, LDA, X, LDX, ONE, B, LDB );

      // Compute the maximum over the number of right hand sides of
      //    norm(B - A*X) / ( max(m,n) * norm(A) * norm(X) * EPS ) .

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 10
         BNORM = dasum( N1, B( 1, J ), 1 );
         XNORM = dasum( N2, X( 1, J ), 1 );
         if ( ANORM == ZERO && BNORM == ZERO ) {
            RESID = ZERO;
         } else if ( ANORM <= ZERO || XNORM <= ZERO ) {
            RESID = ONE / EPS;
         } else {
            RESID = max( RESID, ( ( BNORM / ANORM ) / XNORM ) / ( max( M, N )*EPS ) );
         }
      } // 10

      }
