      REAL sla_gerpvgrw(N, NCOLS, A, LDA, AF, LDAF ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                N, NCOLS, LDA, LDAF;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), AF( LDAF, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, J;
      REAL               AMAX, UMAX, RPVGRW;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      RPVGRW = 1.0;

      for (J = 1; J <= NCOLS; J++) {
         AMAX = 0.0;
         UMAX = 0.0;
         for (I = 1; I <= N; I++) {
            AMAX = max( ( A( I, J ) ).abs(), AMAX );
         }
         for (I = 1; I <= J; I++) {
            UMAX = max( ( AF( I, J ) ).abs(), UMAX );
         }
         if ( UMAX /= 0.0 ) {
            RPVGRW = min( AMAX / UMAX, RPVGRW );
         }
      }
      SLA_GERPVGRW = RPVGRW;
      }
