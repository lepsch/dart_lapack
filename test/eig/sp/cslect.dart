      bool cslect(final int Z,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      Complex            Z;
      // ..

      double               ZERO;
      const              ZERO = 0.0 ;
      int                I;
      double               RMIN, X;
      // ..
      // .. Scalars in Common ..
      int                SELDIM, SELOPT;
      // ..
      // .. Arrays in Common ..
      bool               SELVAL( 20 );
      double               SELWI( 20 ), SELWR( 20 );
      // ..
      // .. Common blocks ..
      // COMMON / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, REAL

      if ( SELOPT == 0 ) {
         CSLECT = ( double( Z ) < ZERO );
      } else {
         RMIN = ABS( Z-CMPLX( SELWR( 1 ), SELWI( 1 ) ) );
         CSLECT = SELVAL( 1 );
         for (I = 2; I <= SELDIM; I++) { // 10
            X = ABS( Z-CMPLX( SELWR( I ), SELWI( I ) ) );
            if ( X <= RMIN ) {
               RMIN = X;
               CSLECT = SELVAL( I );
            }
         } // 10
      }
      }
