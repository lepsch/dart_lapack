      bool zlctes(Z, D ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      Complex         D, Z;
      // ..

// =====================================================================

      // .. Parameters ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      double             ZMAX;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX, SIGN
      // ..
      // .. Executable Statements ..

      if ( D == CZERO ) {
         ZLCTES = ( DBLE( Z ) < ZERO );
      } else {
         if ( DBLE( Z ) == ZERO || DBLE( D ) == ZERO ) {
            ZLCTES = ( SIGN( ONE, DIMAG( Z ) ) != SIGN( ONE, DIMAG( D ) ) );
         } else if ( DIMAG( Z ) == ZERO || DIMAG( D ) == ZERO ) {
            ZLCTES = ( SIGN( ONE, DBLE( Z ) ) != SIGN( ONE, DBLE( D ) ) );
         } else {
            ZMAX = max( ( DBLE( Z ) ).abs(), ( DIMAG( Z ) ) ).abs();
            ZLCTES = ( ( DBLE( Z ) / ZMAX )*DBLE( D )+ ( DIMAG( Z ) / ZMAX )*DIMAG( D ) < ZERO );
         }
      }

      return;
      }
