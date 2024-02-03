      bool             FUNCTION ZLCTES( Z, D );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX*16         D, Z
      // ..

*  =====================================================================

      // .. Parameters ..

      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         CZERO
      const              CZERO = ( 0.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      double             ZMAX;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX, SIGN
      // ..
      // .. Executable Statements ..

      if ( D == CZERO ) {
         ZLCTES = ( DBLE( Z ).LT.ZERO )
      } else {
         if ( DBLE( Z ) == ZERO .OR. DBLE( D ) == ZERO ) {
            ZLCTES = ( SIGN( ONE, DIMAG( Z ) ).NE. SIGN( ONE, DIMAG( D ) ) )
         } else if ( DIMAG( Z ) == ZERO .OR. DIMAG( D ) == ZERO ) {
            ZLCTES = ( SIGN( ONE, DBLE( Z ) ).NE. SIGN( ONE, DBLE( D ) ) )
         } else {
            ZMAX = MAX( ABS( DBLE( Z ) ), ABS( DIMAG( Z ) ) )
            ZLCTES = ( ( DBLE( Z ) / ZMAX )*DBLE( D )+ ( DIMAG( Z ) / ZMAX )*DIMAG( D ).LT.ZERO )
         }
      }

      RETURN

      // End of ZLCTES

      }
