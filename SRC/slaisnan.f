      bool    FUNCTION SLAISNAN( SIN1, SIN2 );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL, INTENT(IN) :: SIN1, SIN2
      // ..

*  =====================================================================

*  .. Executable Statements ..
      SLAISNAN = (SIN1.NE.SIN2)
      RETURN
      END
