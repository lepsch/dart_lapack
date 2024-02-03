      REAL sla_gbrpvgrw(N, KL, KU, NCOLS, AB, LDAB, AFB, LDAFB ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                N, KL, KU, NCOLS, LDAB, LDAFB;
      // ..
      // .. Array Arguments ..
      REAL               AB( LDAB, * ), AFB( LDAFB, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, J, KD;
      REAL               AMAX, UMAX, RPVGRW;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      RPVGRW = 1.0;

      KD = KU + 1;
      for (J = 1; J <= NCOLS; J++) {
         AMAX = 0.0;
         UMAX = 0.0;
         DO I = max( J-KU, 1 ), min( J+KL, N );
            AMAX = max( ABS( AB( KD+I-J, J)), AMAX );
         }
         DO I = max( J-KU, 1 ), J;
            UMAX = max( ABS( AFB( KD+I-J, J ) ), UMAX );
         }
         if ( UMAX /= 0.0 ) {
            RPVGRW = min( AMAX / UMAX, RPVGRW );
         }
      }
      SLA_GBRPVGRW = RPVGRW;
      }
