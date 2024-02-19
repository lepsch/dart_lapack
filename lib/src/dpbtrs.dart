import 'dart:math';

import 'package:lapack/src/blas/dtbsv.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dpbtrs(
  final String UPLO,
  final int N,
  final int KD,
  final int NRHS,
  final Matrix<double> AB_,
  final int LDAB,
  final Matrix<double> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.dim(LDAB);
  final B = B_.dim(LDB);
  bool UPPER;
  int J;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (KD < 0) {
    INFO.value = -3;
  } else if (NRHS < 0) {
    INFO.value = -4;
  } else if (LDAB < KD + 1) {
    INFO.value = -6;
  } else if (LDB < max(1, N)) {
    INFO.value = -8;
  }
  if (INFO.value != 0) {
    xerbla('DPBTRS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0) return;

  if (UPPER) {
    // Solve A*X = B where A = U**T *U.

    for (J = 1; J <= NRHS; J++) {
      // Solve U**T *X = B, overwriting B with X.

      dtbsv('Upper', 'Transpose', 'Non-unit', N, KD, AB, LDAB,
          B(1, J).asArray(), 1);

      // Solve U*X = B, overwriting B with X.

      dtbsv('Upper', 'No transpose', 'Non-unit', N, KD, AB, LDAB,
          B(1, J).asArray(), 1);
    }
  } else {
    // Solve A*X = B where A = L*L**T.

    for (J = 1; J <= NRHS; J++) {
      // Solve L*X = B, overwriting B with X.

      dtbsv('Lower', 'No transpose', 'Non-unit', N, KD, AB, LDAB,
          B(1, J).asArray(), 1);

      // Solve L**T *X = B, overwriting B with X.

      dtbsv('Lower', 'Transpose', 'Non-unit', N, KD, AB, LDAB,
          B(1, J).asArray(), 1);
    }
  }
}
