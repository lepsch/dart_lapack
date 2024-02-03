      bool             FUNCTION CSLECT( Z );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX            Z
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0E0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      REAL               RMIN, X
      // ..
      // .. Scalars in Common ..
      int                SELDIM, SELOPT;
      // ..
      // .. Arrays in Common ..
      bool               SELVAL( 20 );
      REAL               SELWI( 20 ), SELWR( 20 )
      // ..
      // .. Common blocks ..
      COMMON             / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, REAL
      // ..
      // .. Executable Statements ..

      if ( SELOPT.EQ.0 ) {
         CSLECT = ( REAL( Z ).LT.ZERO )
      } else {
         RMIN = ABS( Z-CMPLX( SELWR( 1 ), SELWI( 1 ) ) )
         CSLECT = SELVAL( 1 )
         for (I = 2; I <= SELDIM; I++) { // 10
            X = ABS( Z-CMPLX( SELWR( I ), SELWI( I ) ) )
            if ( X.LE.RMIN ) {
               RMIN = X
               CSLECT = SELVAL( I )
            }
         } // 10
      }
      RETURN

      // End of CSLECT

      }
