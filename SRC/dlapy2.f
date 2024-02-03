      double           FUNCTION DLAPY2( X, Y );
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      double             X, Y;
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO;
      PARAMETER          ( ZERO = 0.0D0 )
      double             ONE;
      PARAMETER          ( ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      double             W, XABS, YABS, Z, HUGEVAL;
      bool               X_IS_NAN, Y_IS_NAN;
*     ..
*     .. External Functions ..
      bool               DISNAN;
      EXTERNAL           DISNAN
*     ..
*     .. External Subroutines ..
      double             DLAMCH;
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      X_IS_NAN = DISNAN( X )
      Y_IS_NAN = DISNAN( Y )
      IF ( X_IS_NAN ) DLAPY2 = X
      IF ( Y_IS_NAN ) DLAPY2 = Y
      HUGEVAL = DLAMCH( 'Overflow' )
*
      IF ( .NOT.( X_IS_NAN.OR.Y_IS_NAN ) ) THEN
         XABS = ABS( X )
         YABS = ABS( Y )
         W = MAX( XABS, YABS )
         Z = MIN( XABS, YABS )
         IF( Z.EQ.ZERO .OR. W.GT.HUGEVAL ) THEN
            DLAPY2 = W
         ELSE
            DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
         END IF
      END IF
      RETURN
*
*     End of DLAPY2
*
      END
