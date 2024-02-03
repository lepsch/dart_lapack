      double dcabs1(Z) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      Complex Z;
      // ..
      // ..
// =====================================================================

      // .. Intrinsic Functions ..
      // INTRINSIC ABS,DBLE,DIMAG

      DCABS1 = (DBLE(Z)).abs() + (DIMAG(Z)).abs();
      return;
      }
