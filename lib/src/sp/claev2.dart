void claev2(A, B, C, RT1, RT2, CS1, SN1) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  REAL CS1, RT1, RT2;
  Complex A, B, C, SN1;
  // ..

// =====================================================================

  // .. Parameters ..
  REAL ZERO;
  const ZERO = 0.0;
  REAL ONE;
  const ONE = 1.0;
  // ..
  // .. Local Scalars ..
  REAL T;
  Complex W;
  // ..
  // .. External Subroutines ..
  // EXTERNAL SLAEV2
  // ..
  // .. Intrinsic Functions ..
  // INTRINSIC ABS, CONJG, REAL
  // ..
  // .. Executable Statements ..

  if ((B).abs() == ZERO) {
    W = ONE;
  } else {
    W = CONJG(B) / (B).abs();
  }
  slaev2(REAL(A), (B).abs(), REAL(C), RT1, RT2, CS1, T);
  SN1 = W * T;
  return;
}