String chla_transtype(int TRANS) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const BLAS_NO_TRANS = 111, BLAS_TRANS = 112, BLAS_CONJ_TRANS = 113;

  if (TRANS == BLAS_NO_TRANS) {
    return 'N';
  } else if (TRANS == BLAS_TRANS) {
    return 'T';
  } else if (TRANS == BLAS_CONJ_TRANS) {
    return 'C';
  } else {
    return 'X';
  }
}
