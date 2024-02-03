      REAL FUNCTION SLA_GERPVGRW( N, NCOLS, A, LDA, AF, LDAF )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                N, NCOLS, LDA, LDAF;
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), AF( LDAF, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      int                I, J;
      REAL               AMAX, UMAX, RPVGRW
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      RPVGRW = 1.0

      DO J = 1, NCOLS
         AMAX = 0.0
         UMAX = 0.0
         DO I = 1, N
            AMAX = MAX( ABS( A( I, J ) ), AMAX )
         END DO
         DO I = 1, J
            UMAX = MAX( ABS( AF( I, J ) ), UMAX )
         END DO
         IF ( UMAX /= 0.0 ) THEN
            RPVGRW = MIN( AMAX / UMAX, RPVGRW )
         END IF
      END DO
      SLA_GERPVGRW = RPVGRW
*
*     End of SLA_GERPVGRW
*
      END
