import 'dart:math';

import 'package:lapack/src/blas/dtrsm.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dpotrs(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  const ONE = 1.0;
  bool UPPER;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDB < max(1, N)) {
    INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('DPOTRS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0) return;

  if (UPPER) {
    // Solve A*X = B where A = U**T *U.

    // Solve U**T *X = B, overwriting B with X.

    dtrsm(
        'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB);

    // Solve U*X = B, overwriting B with X.

    dtrsm('Left', 'Upper', 'No transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B,
        LDB);
  } else {
    // Solve A*X = B where A = L*L**T.

    // Solve L*X = B, overwriting B with X.

    dtrsm('Left', 'Lower', 'No transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B,
        LDB);

    // Solve L**T *X = B, overwriting B with X.

    dtrsm(
        'Left', 'Lower', 'Transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB);
  }
}
