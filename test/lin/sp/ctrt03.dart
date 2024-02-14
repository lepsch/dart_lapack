      void ctrt03(final int UPLO, final int TRANS, final int DIAG, final int N, final int NRHS, final Matrix<double> A_, final int LDA, final int SCALE, final int CNORM, final int TSCAL, final Matrix<double> X_, final int LDX, final Matrix<double> B_, final int LDB, final Array<double> _WORK_, final int RESID,) {
  final A = A_.dim();
  final X = X_.dim();
  final B = B_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, TRANS, UPLO;
      int                LDA, LDB, LDX, N, NRHS;
      double               RESID, SCALE, TSCAL;
      double               CNORM( * );
      Complex            A( LDA, * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                IX, J;
      double               EPS, ERR, SMLNUM, TNORM, XNORM, XSCAL;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ICAMAX;
      //- REAL               SLAMCH;
      // EXTERNAL lsame, ICAMAX, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CSSCAL, CTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, MAX, REAL

      // Quick exit if N = 0

      if ( N <= 0 || NRHS <= 0 ) {
         RESID = ZERO;
         return;
      }
      EPS = SLAMCH( 'Epsilon' );
      SMLNUM = SLAMCH( 'Safe minimum' );

      // Compute the norm of the triangular matrix A using the column
      // norms already computed by CLATRS.

      TNORM = ZERO;
      if ( lsame( DIAG, 'N' ) ) {
         for (J = 1; J <= N; J++) { // 10
            TNORM = max( TNORM, TSCAL*( A( J, J ) ).abs()+CNORM( J ) );
         } // 10
      } else {
         for (J = 1; J <= N; J++) { // 20
            TNORM = max( TNORM, TSCAL+CNORM( J ) );
         } // 20
      }

      // Compute the maximum over the number of right hand sides of
      //    norm(op(A)*x - s*b) / ( norm(op(A)) * norm(x) * EPS ).

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 30
         ccopy(N, X( 1, J ), 1, WORK, 1 );
         IX = ICAMAX( N, WORK, 1 );
         XNORM = max( ONE, ( X( IX, J ) ).abs() );
         XSCAL = ( ONE / XNORM ) / REAL( N );
         csscal(N, XSCAL, WORK, 1 );
         ctrmv(UPLO, TRANS, DIAG, N, A, LDA, WORK, 1 );
         caxpy(N, CMPLX( -SCALE*XSCAL ), B( 1, J ), 1, WORK, 1 );
         IX = ICAMAX( N, WORK, 1 );
         ERR = TSCAL*( WORK( IX ) ).abs();
         IX = ICAMAX( N, X( 1, J ), 1 );
         XNORM = ( X( IX, J ) ).abs();
         if ( ERR*SMLNUM <= XNORM ) {
            if (XNORM > ZERO) ERR = ERR / XNORM;
         } else {
            if (ERR > ZERO) ERR = ONE / EPS;
         }
         if ( ERR*SMLNUM <= TNORM ) {
            if (TNORM > ZERO) ERR = ERR / TNORM;
         } else {
            if (ERR > ZERO) ERR = ONE / EPS;
         }
         RESID = max( RESID, ERR );
      } // 30

      }
