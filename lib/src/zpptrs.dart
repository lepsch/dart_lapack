import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/ztpsv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zpptrs(
  final String UPLO,
  final int N,
  final int NRHS,
  final Array<Complex> AP_,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.dim();
  final B = B_.dim(LDB);
  bool UPPER;
  int I;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDB < max(1, N)) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('ZPPTRS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0) return;

  if (UPPER) {
    // Solve A*X = B where A = U**H * U.

    for (I = 1; I <= NRHS; I++) {
      // 10

      // Solve U**H *X = B, overwriting B with X.

      ztpsv('Upper', 'Conjugate transpose', 'Non-unit', N, AP,
          B(1, I).asArray(), 1);

      // Solve U*X = B, overwriting B with X.

      ztpsv('Upper', 'No transpose', 'Non-unit', N, AP, B(1, I).asArray(), 1);
    } // 10
  } else {
    // Solve A*X = B where A = L * L**H.

    for (I = 1; I <= NRHS; I++) {
      // 20

      // Solve L*Y = B, overwriting B with X.

      ztpsv('Lower', 'No transpose', 'Non-unit', N, AP, B(1, I).asArray(), 1);

      // Solve L**H *X = Y, overwriting B with X.

      ztpsv('Lower', 'Conjugate transpose', 'Non-unit', N, AP,
          B(1, I).asArray(), 1);
    } // 20
  }
}
