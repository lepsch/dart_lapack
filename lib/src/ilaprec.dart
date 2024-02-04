int ilaprec(PREC) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  String PREC;
  // ..

// =====================================================================

  // .. Parameters ..
  int BLAS_PREC_SINGLE, BLAS_PREC_DOUBLE, BLAS_PREC_INDIGENOUS, BLAS_PREC_EXTRA;
  const BLAS_PREC_SINGLE = 211,
      BLAS_PREC_DOUBLE = 212,
      BLAS_PREC_INDIGENOUS = 213,
      BLAS_PREC_EXTRA = 214;
  // ..
  // .. External Functions ..
  //- bool               lsame;
  // EXTERNAL lsame
  // ..
  // .. Executable Statements ..
  if (lsame(PREC, 'S')) {
    ILAPREC = BLAS_PREC_SINGLE;
  } else if (lsame(PREC, 'D')) {
    ILAPREC = BLAS_PREC_DOUBLE;
  } else if (lsame(PREC, 'I')) {
    ILAPREC = BLAS_PREC_INDIGENOUS;
  } else if (lsame(PREC, 'X') || lsame(PREC, 'E')) {
    ILAPREC = BLAS_PREC_EXTRA;
  } else {
    ILAPREC = -1;
  }
  return;
}
