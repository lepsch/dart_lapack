      double           FUNCTION ZLA_GBRPVGRW( N, KL, KU, NCOLS, AB, LDAB, AFB, LDAFB );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                N, KL, KU, NCOLS, LDAB, LDAFB;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         AB( LDAB, * ), AFB( LDAFB, * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, J, KD;
      double             AMAX, UMAX, RPVGRW;
      COMPLEX*16         ZDUM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL, DIMAG
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      RPVGRW = 1.0D+0

      KD = KU + 1
      DO J = 1, NCOLS
         AMAX = 0.0D+0
         UMAX = 0.0D+0
         DO I = MAX( J-KU, 1 ), MIN( J+KL, N )
            AMAX = MAX( CABS1( AB( KD+I-J, J ) ), AMAX )
         END DO
         DO I = MAX( J-KU, 1 ), J
            UMAX = MAX( CABS1( AFB( KD+I-J, J ) ), UMAX )
         END DO
         if ( UMAX /= 0.0D+0 ) {
            RPVGRW = MIN( AMAX / UMAX, RPVGRW )
         }
      END DO
      ZLA_GBRPVGRW = RPVGRW

      // End of ZLA_GBRPVGRW

      }
