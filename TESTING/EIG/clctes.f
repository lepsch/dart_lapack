      bool             FUNCTION CLCTES( Z, D );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX            D, Z
      // ..

*  =====================================================================

      // .. Parameters ..

      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CZERO
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
      // ..
      // .. Local Scalars ..
      REAL               ZMAX
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, REAL, SIGN
      // ..
      // .. Executable Statements ..

      IF( D.EQ.CZERO ) THEN
         CLCTES = ( REAL( Z ).LT.ZERO )
      ELSE
         IF( REAL( Z ).EQ.ZERO .OR. REAL( D ).EQ.ZERO ) THEN
            CLCTES = ( SIGN( ONE, AIMAG( Z ) ).NE. SIGN( ONE, AIMAG( D ) ) )
         ELSE IF( AIMAG( Z ).EQ.ZERO .OR. AIMAG( D ).EQ.ZERO ) THEN
            CLCTES = ( SIGN( ONE, REAL( Z ) ).NE. SIGN( ONE, REAL( D ) ) )
         ELSE
            ZMAX = MAX( ABS( REAL( Z ) ), ABS( AIMAG( Z ) ) )
            CLCTES = ( ( REAL( Z ) / ZMAX )*REAL( D )+ ( AIMAG( Z ) / ZMAX )*AIMAG( D ).LT.ZERO )
         END IF
      END IF

      RETURN

      // End of CLCTES

      }
