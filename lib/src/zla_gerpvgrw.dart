      double zla_gerpvgrw(N, NCOLS, A, LDA, AF, LDAF ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                N, NCOLS, LDA, LDAF;
      // ..
      // .. Array Arguments ..
      Complex         A( LDA, * ), AF( LDAF, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, J;
      double             AMAX, UMAX, RPVGRW;
      Complex         ZDUM;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, ABS, REAL, DIMAG
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function Definitions ..
      CABS1[ZDUM] = ( ZDUM.toDouble() ).abs() + ( DIMAG( ZDUM ) ).abs();
      // ..
      // .. Executable Statements ..

      RPVGRW = 1.0;

      for (J = 1; J <= NCOLS; J++) {
         AMAX = 0.0;
         UMAX = 0.0;
         for (I = 1; I <= N; I++) {
            AMAX = max( CABS1( A( I, J ) ), AMAX );
         }
         for (I = 1; I <= J; I++) {
            UMAX = max( CABS1( AF( I, J ) ), UMAX );
         }
         if ( UMAX /= 0.0 ) {
            RPVGRW = min( AMAX / UMAX, RPVGRW );
         }
      }
      ZLA_GERPVGRW = RPVGRW;
      }
