double droundup_lwork(LWORK) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  int LWORK;
  // ..

// =====================================================================
  // ..
  // .. Intrinsic Functions ..
  // INTRINSIC EPSILON, DBLE, INT
  // ..
  // .. Executable Statements ..
  // ..
  DROUNDUP_LWORK = LWORK.toDouble();

  if (INT(DROUNDUP_LWORK) < LWORK) {
    // Force round up of LWORK
    DROUNDUP_LWORK = DROUNDUP_LWORK * (1.0 + EPSILON(0.0));
  }
}
