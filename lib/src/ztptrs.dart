import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/ztpsv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void ztptrs(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int N,
  final int NRHS,
  final Array<Complex> AP_,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> INFO,
) {
  final AP = AP_.having();
  final B = B_.having(ld: LDB);

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  bool NOUNIT, UPPER;
  int J, JC;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  NOUNIT = lsame(DIAG, 'N');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (!lsame(TRANS, 'N') && !lsame(TRANS, 'T') && !lsame(TRANS, 'C')) {
    INFO.value = -2;
  } else if (!NOUNIT && !lsame(DIAG, 'U')) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (NRHS < 0) {
    INFO.value = -5;
  } else if (LDB < max(1, N)) {
    INFO.value = -8;
  }
  if (INFO.value != 0) {
    xerbla('ZTPTRS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Check for singularity.

  if (NOUNIT) {
    if (UPPER) {
      JC = 1;
      for (INFO.value = 1; INFO.value <= N; INFO.value++) {
        if (AP[JC + INFO.value - 1] == Complex.zero) return;
        JC += INFO.value;
      }
    } else {
      JC = 1;
      for (INFO.value = 1; INFO.value <= N; INFO.value++) {
        if (AP[JC] == Complex.zero) return;
        JC += N - INFO.value + 1;
      }
    }
  }
  INFO.value = 0;

  // Solve  A * x = b,  A**T * x = b,  or  A**H * x = b.

  for (J = 1; J <= NRHS; J++) {
    ztpsv(UPLO, TRANS, DIAG, N, AP, B(1, J).asArray(), 1);
  }
}
