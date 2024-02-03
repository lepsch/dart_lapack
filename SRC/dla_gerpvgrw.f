      double           FUNCTION DLA_GERPVGRW( N, NCOLS, A, LDA, AF, LDAF );
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                N, NCOLS, LDA, LDAF;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), AF( LDAF, * );
      // ..
*
*  =====================================================================
*
      // .. Local Scalars ..
      int                I, J;
      double             AMAX, UMAX, RPVGRW;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..
*
      RPVGRW = 1.0D+0

      DO J = 1, NCOLS
         AMAX = 0.0D+0
         UMAX = 0.0D+0
         DO I = 1, N
            AMAX = MAX( ABS( A( I, J ) ), AMAX )
         END DO
         DO I = 1, J
            UMAX = MAX( ABS( AF( I, J ) ), UMAX )
         END DO
         IF ( UMAX /= 0.0D+0 ) THEN
            RPVGRW = MIN( AMAX / UMAX, RPVGRW )
         END IF
      END DO
      DLA_GERPVGRW = RPVGRW
*
      // End of DLA_GERPVGRW
*
      END
