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
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            CZERO
      const              CZERO = ( 0.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      REAL               ZMAX
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, REAL, SIGN
      // ..
      // .. Executable Statements ..

      if ( D == CZERO ) {
         CLCTES = ( REAL( Z ).LT.ZERO )
      } else {
         if ( REAL( Z ) == ZERO .OR. REAL( D ) == ZERO ) {
            CLCTES = ( SIGN( ONE, AIMAG( Z ) ).NE. SIGN( ONE, AIMAG( D ) ) )
         } else if ( AIMAG( Z ) == ZERO .OR. AIMAG( D ) == ZERO ) {
            CLCTES = ( SIGN( ONE, REAL( Z ) ).NE. SIGN( ONE, REAL( D ) ) )
         } else {
            ZMAX = MAX( ABS( REAL( Z ) ), ABS( AIMAG( Z ) ) )
            CLCTES = ( ( REAL( Z ) / ZMAX )*REAL( D )+ ( AIMAG( Z ) / ZMAX )*AIMAG( D ).LT.ZERO )
         }
      }

      RETURN

      // End of CLCTES

      }
