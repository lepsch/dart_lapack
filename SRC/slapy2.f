      REAL             FUNCTION SLAPY2( X, Y )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL               X, Y
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      REAL               W, XABS, YABS, Z, HUGEVAL
      bool               X_IS_NAN, Y_IS_NAN;
*     ..
*     .. External Functions ..
      bool               SISNAN;
      EXTERNAL           SISNAN
*     ..
*     .. External Subroutines ..
      REAL               SLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      X_IS_NAN = SISNAN( X )
      Y_IS_NAN = SISNAN( Y )
      IF ( X_IS_NAN ) SLAPY2 = X
      IF ( Y_IS_NAN ) SLAPY2 = Y
      HUGEVAL = SLAMCH( 'Overflow' )
*
      IF ( .NOT.( X_IS_NAN.OR.Y_IS_NAN ) ) THEN
         XABS = ABS( X )
         YABS = ABS( Y )
         W = MAX( XABS, YABS )
         Z = MIN( XABS, YABS )
         IF( Z.EQ.ZERO .OR. W.GT.HUGEVAL ) THEN
            SLAPY2 = W
         ELSE
            SLAPY2 = W*SQRT( ONE+( Z / W )**2 )
         END IF
      END IF
      RETURN
*
*     End of SLAPY2
*
      END
