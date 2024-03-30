import 'dart:math';

import 'package:lapack/src/blas/dtbsv.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dtbtrs(
  final String UPLO,
  final String TRANS,
  final String DIAG,
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
  final AB = AB_.having(ld: LDAB);
  final B = B_.having(ld: LDB);
  const ZERO = 0.0;
  bool NOUNIT, UPPER;
  int J;

  // Test the input parameters.

  INFO.value = 0;
  NOUNIT = lsame(DIAG, 'N');
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (!lsame(TRANS, 'N') && !lsame(TRANS, 'T') && !lsame(TRANS, 'C')) {
    INFO.value = -2;
  } else if (!NOUNIT && !lsame(DIAG, 'U')) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (KD < 0) {
    INFO.value = -5;
  } else if (NRHS < 0) {
    INFO.value = -6;
  } else if (LDAB < KD + 1) {
    INFO.value = -8;
  } else if (LDB < max(1, N)) {
    INFO.value = -10;
  }
  if (INFO.value != 0) {
    xerbla('DTBTRS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Check for singularity.

  if (NOUNIT) {
    if (UPPER) {
      for (INFO.value = 1; INFO.value <= N; INFO.value++) {
        if (AB[KD + 1][INFO.value] == ZERO) return;
      }
    } else {
      for (INFO.value = 1; INFO.value <= N; INFO.value++) {
        if (AB[1][INFO.value] == ZERO) return;
      }
    }
  }
  INFO.value = 0;

  // Solve A * X = B  or  A**T * X = B.

  for (J = 1; J <= NRHS; J++) {
    dtbsv(UPLO, TRANS, DIAG, N, KD, AB, LDAB, B(1, J).asArray(), 1);
  }
}
