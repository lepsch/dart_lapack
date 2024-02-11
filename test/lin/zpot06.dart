      void zpot06(final int UPLO, final int N, final int NRHS, final Matrix<double> A, final int LDA, final Matrix<double> X, final int LDX, final Matrix<double> B, final int LDB, final Array<double> RWORK, final int RESID,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                LDA, LDB, LDX, N, NRHS;
      double             RESID;
      double             RWORK( * );
      Complex         A( LDA, * ), B( LDB, * ), X( LDX, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CONE, NEGCONE;
      const              CONE = ( 1.0, 0.0 ) ;
      const              NEGCONE = ( -1.0, 0.0 ) ;
      int                IFAIL, J;
      double             ANORM, BNORM, EPS, XNORM;
      Complex         ZDUM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                IZAMAX;
      //- double             DLAMCH, ZLANSY;
      // EXTERNAL lsame, IZAMAX, DLAMCH, ZLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZHEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[ZDUM] = ( ZDUM.toDouble() ).abs() + ( DIMAG( ZDUM ) ).abs();
      // ..

      // Quick exit if N = 0 or NRHS = 0

      if ( N <= 0 || NRHS == 0 ) {
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = dlamch( 'Epsilon' );
      ANORM = ZLANSY( 'I', UPLO, N, A, LDA, RWORK );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Compute  B - A*X  and store in B.
      IFAIL=0;

      zhemm('Left', UPLO, N, NRHS, NEGCONE, A, LDA, X, LDX, CONE, B, LDB );

      // Compute the maximum over the number of right hand sides of
      //    norm(B - A*X) / ( norm(A) * norm(X) * EPS ) .

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 10
         BNORM = CABS1(B(IZAMAX( N, B( 1, J ), 1 ),J));
         XNORM = CABS1(X(IZAMAX( N, X( 1, J ), 1 ),J));
         if ( XNORM <= ZERO ) {
            RESID = ONE / EPS;
         } else {
            RESID = max( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS );
         }
      } // 10

      }
