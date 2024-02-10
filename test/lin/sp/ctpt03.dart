      void ctpt03(final int UPLO, final int TRANS, final int DIAG, final int N, final int NRHS, final int AP, final int SCALE, final int CNORM, final int TSCAL, final Matrix<double> X, final int LDX, final Matrix<double> B, final int LDB, final Array<double> _WORK, final int RESID) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, TRANS, UPLO;
      int                LDB, LDX, N, NRHS;
      double               RESID, SCALE, TSCAL;
      double               CNORM( * );
      Complex            AP( * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                IX, J, JJ;
      double               EPS, ERR, SMLNUM, TNORM, XNORM, XSCAL;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ICAMAX;
      //- REAL               SLAMCH;
      // EXTERNAL lsame, ICAMAX, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CSSCAL, CTPMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, MAX, REAL

      // Quick exit if N = 0.

      if ( N <= 0 || NRHS <= 0 ) {
         RESID = ZERO;
         return;
      }
      EPS = SLAMCH( 'Epsilon' );
      SMLNUM = SLAMCH( 'Safe minimum' );

      // Compute the norm of the triangular matrix A using the column
      // norms already computed by CLATPS.

      TNORM = 0.;
      if ( lsame( DIAG, 'N' ) ) {
         if ( lsame( UPLO, 'U' ) ) {
            JJ = 1;
            for (J = 1; J <= N; J++) { // 10
               TNORM = max( TNORM, TSCAL*( AP( JJ ) ).abs()+CNORM( J ) );
               JJ = JJ + J + 1;
            } // 10
         } else {
            JJ = 1;
            for (J = 1; J <= N; J++) { // 20
               TNORM = max( TNORM, TSCAL*( AP( JJ ) ).abs()+CNORM( J ) );
               JJ = JJ + N - J + 1;
            } // 20
         }
      } else {
         for (J = 1; J <= N; J++) { // 30
            TNORM = max( TNORM, TSCAL+CNORM( J ) );
         } // 30
      }

      // Compute the maximum over the number of right hand sides of
      //    norm(op(A)*x - s*b) / ( norm(A) * norm(x) * EPS ).

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 40
         ccopy(N, X( 1, J ), 1, WORK, 1 );
         IX = ICAMAX( N, WORK, 1 );
         XNORM = max( ONE, ( X( IX, J ) ).abs() );
         XSCAL = ( ONE / XNORM ) / REAL( N );
         csscal(N, XSCAL, WORK, 1 );
         ctpmv(UPLO, TRANS, DIAG, N, AP, WORK, 1 );
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
      } // 40

      }
