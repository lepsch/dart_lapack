      double sla_gbrpvgrw(N, KL, KU, NCOLS, final Matrix<double> AB, final int LDAB, AFB, LDAFB ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                N, KL, KU, NCOLS, LDAB, LDAFB;
      double               AB( LDAB, * ), AFB( LDAFB, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, J, KD;
      double               AMAX, UMAX, RPVGRW;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN

      RPVGRW = 1.0;

      KD = KU + 1;
      for (J = 1; J <= NCOLS; J++) {
         AMAX = 0.0;
         UMAX = 0.0;
         for (I = max( J-KU, 1 ); I <= min( J+KL, N ); I++) {
            AMAX = max( ( AB( KD+I-J, J)).abs(), AMAX );
         }
         for (I = max( J-KU, 1 ); I <= J; I++) {
            UMAX = max( ( AFB( KD+I-J, J ) ).abs(), UMAX );
         }
         if ( UMAX /= 0.0 ) {
            RPVGRW = min( AMAX / UMAX, RPVGRW );
         }
      }
      SLA_GBRPVGRW = RPVGRW;
      }
