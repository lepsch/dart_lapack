      REAL FUNCTION SLA_GBRPVGRW( N, KL, KU, NCOLS, AB, LDAB, AFB, LDAFB )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                N, KL, KU, NCOLS, LDAB, LDAFB;
      // ..
      // .. Array Arguments ..
      REAL               AB( LDAB, * ), AFB( LDAFB, * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, J, KD;
      REAL               AMAX, UMAX, RPVGRW
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      RPVGRW = 1.0

      KD = KU + 1
      DO J = 1, NCOLS
         AMAX = 0.0
         UMAX = 0.0
         DO I = MAX( J-KU, 1 ), MIN( J+KL, N )
            AMAX = MAX( ABS( AB( KD+I-J, J)), AMAX )
         END DO
         DO I = MAX( J-KU, 1 ), J
            UMAX = MAX( ABS( AFB( KD+I-J, J ) ), UMAX )
         END DO
         if ( UMAX /= 0.0 ) {
            RPVGRW = MIN( AMAX / UMAX, RPVGRW )
         }
      END DO
      SLA_GBRPVGRW = RPVGRW

      // End of SLA_GBRPVGRW

      }
