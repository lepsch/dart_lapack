      bool    FUNCTION DISNAN( DIN );
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      double          , INTENT(IN) :: DIN;
      // ..
*
*  =====================================================================
*
*  .. External Functions ..
      bool    DLAISNAN;
      // EXTERNAL DLAISNAN
*  ..
*  .. Executable Statements ..
      DISNAN = DLAISNAN(DIN,DIN)
      RETURN
      END
