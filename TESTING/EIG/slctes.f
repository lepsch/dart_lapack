      bool             FUNCTION SLCTES( ZR, ZI, D );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               D, ZI, ZR
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SIGN
      // ..
      // .. Executable Statements ..

      if ( D == ZERO ) {
         SLCTES = ( ZR.LT.ZERO )
      } else {
         SLCTES = ( SIGN( ONE, ZR ) != SIGN( ONE, D ) )
      }

      RETURN

      // End of SLCTES

      }
