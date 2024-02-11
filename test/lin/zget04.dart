      void zget04(final int N, final int NRHS, final Matrix<double> X, final int LDX, final Matrix<double> XACT, final int LDXACT, final int RCOND, final int RESID,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDX, LDXACT, N, NRHS;
      double             RCOND, RESID;
      Complex         X( LDX, * ), XACT( LDXACT, * );
      // ..

      double             ZERO;
      const              ZERO = 0.0 ;
      int                I, IX, J;
      double             DIFFNM, EPS, XNORM;
      Complex         ZDUM;
      // ..
      // .. External Functions ..
      //- int                IZAMAX;
      //- double             DLAMCH;
      // EXTERNAL IZAMAX, DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[ZDUM] = ( ZDUM.toDouble() ).abs() + ( DIMAG( ZDUM ) ).abs();

      // Quick exit if N = 0 or NRHS = 0.

      if ( N <= 0 || NRHS <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if RCOND is invalid.

      EPS = dlamch( 'Epsilon' );
      if ( RCOND < ZERO ) {
         RESID = 1.0 / EPS;
         return;
      }

      // Compute the maximum of
      //    norm(X - XACT) / ( norm(XACT) * EPS )
      // over all the vectors X and XACT .

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 20
         IX = IZAMAX( N, XACT( 1, J ), 1 );
         XNORM = CABS1( XACT( IX, J ) );
         DIFFNM = ZERO;
         for (I = 1; I <= N; I++) { // 10
            DIFFNM = max( DIFFNM, CABS1( X( I, J )-XACT( I, J ) ) );
         } // 10
         if ( XNORM <= ZERO ) {
            if (DIFFNM > ZERO) RESID = 1.0 / EPS;
         } else {
            RESID = max( RESID, ( DIFFNM / XNORM )*RCOND );
         }
      } // 20
      if (RESID*EPS < 1.0) RESID = RESID / EPS;

      }
