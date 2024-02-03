      void dget04(N, NRHS, X, LDX, XACT, LDXACT, RCOND, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDX, LDXACT, N, NRHS;
      double             RCOND, RESID;
      // ..
      // .. Array Arguments ..
      double             X( LDX, * ), XACT( LDXACT, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I, IX, J;
      double             DIFFNM, EPS, XNORM;
      // ..
      // .. External Functions ..
      int                IDAMAX;
      double             DLAMCH;
      // EXTERNAL IDAMAX, DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0.

      if ( N <= 0 || NRHS <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if RCOND is invalid.

      EPS = DLAMCH( 'Epsilon' );
      if ( RCOND < ZERO ) {
         RESID = 1.0 / EPS;
         return;
      }

      // Compute the maximum of
         // norm(X - XACT) / ( norm(XACT) * EPS )
      // over all the vectors X and XACT .

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 20
         IX = IDAMAX( N, XACT( 1, J ), 1 );
         XNORM = ABS( XACT( IX, J ) );
         DIFFNM = ZERO;
         for (I = 1; I <= N; I++) { // 10
            DIFFNM = max( DIFFNM, ABS( X( I, J )-XACT( I, J ) ) );
         } // 10
         if ( XNORM <= ZERO ) {
            if (DIFFNM > ZERO) RESID = 1.0 / EPS;
         } else {
            RESID = max( RESID, ( DIFFNM / XNORM )*RCOND );
         }
      } // 20
      if (RESID*EPS < 1.0) RESID = RESID / EPS;

      return;

      // End of DGET04

      }
