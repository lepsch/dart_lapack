int ilauplo(UPLO) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  String UPLO;
  // ..

// =====================================================================

  // .. Parameters ..
  int BLAS_UPPER, BLAS_LOWER;
  const BLAS_UPPER = 121, BLAS_LOWER = 122;
  // ..
  // .. External Functions ..
  //- bool               lsame;
  // EXTERNAL lsame
  // ..
  // .. Executable Statements ..
  if (lsame(UPLO, 'U')) {
    ILAUPLO = BLAS_UPPER;
  } else if (lsame(UPLO, 'L')) {
    ILAUPLO = BLAS_LOWER;
  } else {
    ILAUPLO = -1;
  }
  return;
}
