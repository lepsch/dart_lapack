import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dpptrf.dart';
import 'package:lapack/src/dpptrs.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dppsv(
  final String UPLO,
  final int N,
  final int NRHS,
  final Array<double> AP_,
  final Matrix<double> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.dim();
  final B = B_.dim(LDB);

  // Test the input parameters.

  INFO.value = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDB < max(1, N)) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('DPPSV ', -INFO.value);
    return;
  }

  // Compute the Cholesky factorization A = U**T*U or A = L*L**T.

  dpptrf(UPLO, N, AP, INFO);
  if (INFO.value == 0) {
    // Solve the system A*X = B, overwriting B with X.

    dpptrs(UPLO, N, NRHS, AP, B, LDB, INFO);
  }
}
