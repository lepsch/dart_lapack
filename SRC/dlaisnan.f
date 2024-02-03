      bool    FUNCTION DLAISNAN( DIN1, DIN2 );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double          , INTENT(IN) :: DIN1, DIN2;
      // ..

*  =====================================================================

*  .. Executable Statements ..
      DLAISNAN = (DIN1.NE.DIN2)
      RETURN
      END
