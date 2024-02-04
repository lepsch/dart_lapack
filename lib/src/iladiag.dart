int iladiag(DIAG) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  String DIAG;
  // ..

// =====================================================================

  // .. Parameters ..
  int BLAS_NON_UNIT_DIAG, BLAS_UNIT_DIAG;
  const BLAS_NON_UNIT_DIAG = 131, BLAS_UNIT_DIAG = 132;
  // ..
  // .. External Functions ..
  //- bool               lsame;
  // EXTERNAL lsame
  // ..
  // .. Executable Statements ..
  if (lsame(DIAG, 'N')) {
    ILADIAG = BLAS_NON_UNIT_DIAG;
  } else if (lsame(DIAG, 'U')) {
    ILADIAG = BLAS_UNIT_DIAG;
  } else {
    ILADIAG = -1;
  }
  return;
}
