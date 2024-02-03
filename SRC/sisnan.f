      bool    FUNCTION SISNAN( SIN );
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL, INTENT(IN) :: SIN
*     ..
*
*  =====================================================================
*
*  .. External Functions ..
      bool    SLAISNAN;
      EXTERNAL SLAISNAN
*  ..
*  .. Executable Statements ..
      SISNAN = SLAISNAN(SIN,SIN)
      RETURN
      END
