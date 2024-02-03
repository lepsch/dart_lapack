      bool             FUNCTION ZSLECT( Z );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX*16         Z
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D0 ;
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
      COMMON             / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX
      // ..
      // .. Executable Statements ..

      if ( SELOPT.EQ.0 ) {
         ZSLECT = ( DBLE( Z ).LT.ZERO )
      } else {
         RMIN = ABS( Z-DCMPLX( SELWR( 1 ), SELWI( 1 ) ) )
         ZSLECT = SELVAL( 1 )
         DO 10 I = 2, SELDIM
            X = ABS( Z-DCMPLX( SELWR( I ), SELWI( I ) ) )
            if ( X.LE.RMIN ) {
               RMIN = X
               ZSLECT = SELVAL( I )
            }
   10    CONTINUE
      }
      RETURN

      // End of ZSLECT

      }
