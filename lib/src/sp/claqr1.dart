      void claqr1(N, H, LDH, S1, S2, V ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      Complex            S1, S2;
      int                LDH, N;
      // ..
      // .. Array Arguments ..
      Complex            H( LDH, * ), V( * );
      // ..

// ================================================================

      // .. Parameters ..
      Complex            ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      double               RZERO;
      const              RZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      Complex            CDUM, H21S, H31S;
      double               S;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, REAL
      // ..
      // .. Statement Functions ..
      double               CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[CDUM] = ( double( CDUM ) ).abs() + ( AIMAG( CDUM ) ).abs();
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N != 2 && N != 3 ) {
         return;
      }

      if ( N == 2 ) {
         S = CABS1( H( 1, 1 )-S2 ) + CABS1( H( 2, 1 ) );
         if ( S == RZERO ) {
            V[1] = ZERO;
            V[2] = ZERO;
         } else {
            H21S = H( 2, 1 ) / S;
            V[1] = H21S*H( 1, 2 ) + ( H( 1, 1 )-S1 )* ( ( H( 1, 1 )-S2 ) / S );
            V[2] = H21S*( H( 1, 1 )+H( 2, 2 )-S1-S2 );
         }
      } else {
         S = CABS1( H( 1, 1 )-S2 ) + CABS1( H( 2, 1 ) ) + CABS1( H( 3, 1 ) );
         if ( S == ZERO ) {
            V[1] = ZERO;
            V[2] = ZERO;
            V[3] = ZERO;
         } else {
            H21S = H( 2, 1 ) / S;
            H31S = H( 3, 1 ) / S;
            V[1] = ( H( 1, 1 )-S1 )*( ( H( 1, 1 )-S2 ) / S ) + H( 1, 2 )*H21S + H( 1, 3 )*H31S;
            V[2] = H21S*( H( 1, 1 )+H( 2, 2 )-S1-S2 ) + H( 2, 3 )*H31S;
            V[3] = H31S*( H( 1, 1 )+H( 3, 3 )-S1-S2 ) + H21S*H( 3, 2 );
         }
      }
      }
