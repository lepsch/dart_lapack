      bool clctes(Z, D ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX            D, Z;
      // ..

// =====================================================================

      // .. Parameters ..

      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX            CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      REAL               ZMAX;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, REAL, SIGN
      // ..
      // .. Executable Statements ..

      if ( D == CZERO ) {
         CLCTES = ( REAL( Z ) < ZERO );
      } else {
         if ( REAL( Z ) == ZERO || REAL( D ) == ZERO ) {
            CLCTES = ( SIGN( ONE, AIMAG( Z ) ) != SIGN( ONE, AIMAG( D ) ) );
         } else if ( AIMAG( Z ) == ZERO || AIMAG( D ) == ZERO ) {
            CLCTES = ( SIGN( ONE, REAL( Z ) ) != SIGN( ONE, REAL( D ) ) );
         } else {
            ZMAX = max( ( REAL( Z ) ).abs(), ( AIMAG( Z ) ) ).abs();
            CLCTES = ( ( REAL( Z ) / ZMAX )*REAL( D )+ ( AIMAG( Z ) / ZMAX )*AIMAG( D ) < ZERO );
         }
      }

      return;
      }
