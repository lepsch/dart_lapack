      double           FUNCTION DCABS1(Z);

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX*16 Z
      // ..
      // ..
*  =====================================================================

      // .. Intrinsic Functions ..
      // INTRINSIC ABS,DBLE,DIMAG

      DCABS1 = ABS(DBLE(Z)) + ABS(DIMAG(Z))
      RETURN

      // End of DCABS1

      END
