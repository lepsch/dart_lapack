      bool dslect(ZR, ZI ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             ZI, ZR;
      // ..

// =====================================================================

      // .. Arrays in Common ..
      bool               SELVAL( 20 );
      double             SELWI( 20 ), SELWR( 20 );
      // ..
      // .. Scalars in Common ..
      int                SELDIM, SELOPT;
      // ..
      // .. Common blocks ..
      // COMMON / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
      // ..
      // .. Local Scalars ..
      int                I;
      double             RMIN, X;
      // ..
      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. External Functions ..
      //- double             DLAPY2;
      // EXTERNAL DLAPY2
      // ..
      // .. Executable Statements ..

      if ( SELOPT == 0 ) {
         DSLECT = ( ZR < ZERO );
      } else {
         RMIN = DLAPY2( ZR-SELWR( 1 ), ZI-SELWI( 1 ) );
         DSLECT = SELVAL( 1 );
         for (I = 2; I <= SELDIM; I++) { // 10
            X = DLAPY2( ZR-SELWR( I ), ZI-SELWI( I ) );
            if ( X <= RMIN ) {
               RMIN = X;
               DSLECT = SELVAL( I );
            }
         } // 10
      }
      return;
      }
