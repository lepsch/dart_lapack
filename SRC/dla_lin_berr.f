      void dla_lin_berr(N, NZ, NRHS, RES, AYB, BERR ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                N, NZ, NRHS;
      // ..
      // .. Array Arguments ..
      double             AYB( N, NRHS ), BERR( NRHS );
      double             RES( N, NRHS );
      // ..

// =====================================================================

      // .. Local Scalars ..
      double             TMP;
      int                I, J;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. External Functions ..
      // EXTERNAL DLAMCH
      double             DLAMCH;
      double             SAFE1;
      // ..
      // .. Executable Statements ..

      // Adding SAFE1 to the numerator guards against spuriously zero
      // residuals.  A similar safeguard is in the SLA_yyAMV routine used
      // to compute AYB.

      SAFE1 = DLAMCH( 'Safe minimum' );
      SAFE1 = (NZ+1)*SAFE1;

      for (J = 1; J <= NRHS; J++) {
         BERR(J) = 0.0;
         for (I = 1; I <= N; I++) {
            if (AYB(I,J) != 0.0) {
               TMP = (SAFE1+ABS(RES(I,J)))/AYB(I,J);
               BERR(J) = max( BERR(J), TMP );
            }

      // If AYB is exactly 0.0 (and if computed by SLA_yyAMV), then we know
      // the true residual also must be exactly 0.0.

         }
      }
      }
