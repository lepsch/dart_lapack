      double           FUNCTION ZLA_GERPVGRW( N, NCOLS, A, LDA, AF, LDAF );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                N, NCOLS, LDA, LDAF;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), AF( LDAF, * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, J;
      double             AMAX, UMAX, RPVGRW;
      COMPLEX*16         ZDUM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, ABS, REAL, DIMAG
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      RPVGRW = 1.0D+0

      for (J = 1; J <= NCOLS; J++) {
         AMAX = 0.0D+0
         UMAX = 0.0D+0
         for (I = 1; I <= N; I++) {
            AMAX = MAX( CABS1( A( I, J ) ), AMAX )
         END DO
         for (I = 1; I <= J; I++) {
            UMAX = MAX( CABS1( AF( I, J ) ), UMAX )
         END DO
         if ( UMAX /= 0.0D+0 ) {
            RPVGRW = MIN( AMAX / UMAX, RPVGRW )
         }
      END DO
      ZLA_GERPVGRW = RPVGRW

      // End of ZLA_GERPVGRW

      }
