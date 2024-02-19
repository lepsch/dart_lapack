import 'package:lapack/src/blas/lsame.dart';

int ilauplo(final String UPLO) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const BLAS_UPPER = 121, BLAS_LOWER = 122;

  if (lsame(UPLO, 'U')) {
    return BLAS_UPPER;
  } else if (lsame(UPLO, 'L')) {
    return BLAS_LOWER;
  } else {
    return -1;
  }
}
