      void dtrt03(UPLO, TRANS, DIAG, N, NRHS, final Matrix<double> A, final int LDA, SCALE, CNORM, TSCAL, final Matrix<double> X, final int LDX, final Matrix<double> B, final int LDB, WORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, TRANS, UPLO;
      int                LDA, LDB, LDX, N, NRHS;
      double             RESID, SCALE, TSCAL;
      double             A( LDA, * ), B( LDB, * ), CNORM( * ), WORK( * ), X( LDX, * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                IX, J;
      double             BIGNUM, EPS, ERR, SMLNUM, TNORM, XNORM, XSCAL;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                idamax;
      //- double             DLAMCH;
      // EXTERNAL lsame, idamax, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DSCAL, DTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX

      // Quick exit if N = 0

      if ( N <= 0 || NRHS <= 0 ) {
         RESID = ZERO;
         return;
      }
      EPS = dlamch( 'Epsilon' );
      SMLNUM = dlamch( 'Safe minimum' );
      BIGNUM = ONE / SMLNUM;

      // Compute the norm of the triangular matrix A using the column
      // norms already computed by DLATRS.

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
         dcopy(N, X( 1, J ), 1, WORK, 1 );
         IX = idamax( N, WORK, 1 );
         XNORM = max( ONE, ( X( IX, J ) ).abs() );
         XSCAL = ( ONE / XNORM ) / N.toDouble();
         dscal(N, XSCAL, WORK, 1 );
         dtrmv(UPLO, TRANS, DIAG, N, A, LDA, WORK, 1 );
         daxpy(N, -SCALE*XSCAL, B( 1, J ), 1, WORK, 1 );
         IX = idamax( N, WORK, 1 );
         ERR = TSCAL*( WORK( IX ) ).abs();
         IX = idamax( N, X( 1, J ), 1 );
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
