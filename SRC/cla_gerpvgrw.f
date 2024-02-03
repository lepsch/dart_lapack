      REAL FUNCTION CLA_GERPVGRW( N, NCOLS, A, LDA, AF, LDAF )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                N, NCOLS, LDA, LDAF;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), AF( LDAF, * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, J;
      REAL               AMAX, UMAX, RPVGRW
      COMPLEX            ZDUM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, ABS, REAL, AIMAG
      // ..
      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      RPVGRW = 1.0

      DO J = 1, NCOLS
         AMAX = 0.0
         UMAX = 0.0
         DO I = 1, N
            AMAX = MAX( CABS1( A( I, J ) ), AMAX )
         END DO
         DO I = 1, J
            UMAX = MAX( CABS1( AF( I, J ) ), UMAX )
         END DO
         if ( UMAX /= 0.0 ) {
            RPVGRW = MIN( AMAX / UMAX, RPVGRW )
         }
      END DO
      CLA_GERPVGRW = RPVGRW

      // End of CLA_GERPVGRW

      }
