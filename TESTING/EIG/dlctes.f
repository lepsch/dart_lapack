      bool             FUNCTION DLCTES( ZR, ZI, D );
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      double             D, ZI, ZR;
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SIGN
      // ..
      // .. Executable Statements ..
*
      IF( D.EQ.ZERO ) THEN
         DLCTES = ( ZR.LT.ZERO )
      ELSE
         DLCTES = ( SIGN( ONE, ZR ).NE.SIGN( ONE, D ) )
      END IF
*
      RETURN
*
      // End of DLCTES
*
      END
