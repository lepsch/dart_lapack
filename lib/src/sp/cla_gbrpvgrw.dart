      double cla_gbrpvgrw(final int N, final int KL, final int KU, final int NCOLS, final Matrix<double> AB_, final int LDAB, final int AFB, final int LDAFB,) {
  final AB = AB_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                N, KL, KU, NCOLS, LDAB, LDAFB;
      Complex            AB( LDAB, * ), AFB( LDAFB, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, J, KD;
      double               AMAX, UMAX, RPVGRW;
      Complex            ZDUM;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL, AIMAG
      // ..
      // .. Statement Functions ..
      double               CABS1;
      // ..
      // .. Statement Function Definitions ..
      CABS1[ZDUM] = ( double( ZDUM ) ).abs() + ( AIMAG( ZDUM ) ).abs();

      RPVGRW = 1.0;

      KD = KU + 1;
      for (J = 1; J <= NCOLS; J++) {
         AMAX = 0.0;
         UMAX = 0.0;
         for (I = max( J-KU, 1 ); I <= min( J+KL, N ); I++) {
            AMAX = max( CABS1( AB( KD+I-J, J ) ), AMAX );
         }
         for (I = max( J-KU, 1 ); I <= J; I++) {
            UMAX = max( CABS1( AFB( KD+I-J, J ) ), UMAX );
         }
         if ( UMAX /= 0.0 ) {
            RPVGRW = min( AMAX / UMAX, RPVGRW );
         }
      }
      CLA_GBRPVGRW = RPVGRW;
      }
