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

      IF( D.EQ.ZERO ) THEN
         SLCTES = ( ZR.LT.ZERO )
      } else {
         SLCTES = ( SIGN( ONE, ZR ).NE.SIGN( ONE, D ) )
      END IF

      RETURN

      // End of SLCTES

      }
