      bool zslect(Z ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      Complex         Z;
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      double             RMIN, X;
      // ..
      // .. Scalars in Common ..
      int                SELDIM, SELOPT;
      // ..
      // .. Arrays in Common ..
      bool               SELVAL( 20 );
      double             SELWI( 20 ), SELWR( 20 );
      // ..
      // .. Common blocks ..
      // COMMON / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX
      // ..
      // .. Executable Statements ..

      if ( SELOPT == 0 ) {
         ZSLECT = ( Z.toDouble() < ZERO );
      } else {
         RMIN = ABS( Z-DCMPLX( SELWR( 1 ), SELWI( 1 ) ) );
         ZSLECT = SELVAL( 1 );
         for (I = 2; I <= SELDIM; I++) { // 10
            X = ABS( Z-DCMPLX( SELWR( I ), SELWI( I ) ) );
            if ( X <= RMIN ) {
               RMIN = X;
               ZSLECT = SELVAL( I );
            }
         } // 10
      }
      return;
      }
