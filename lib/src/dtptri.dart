import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dtpmv.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dtptri(
  final String UPLO,
  final String DIAG,
  final int N,
  final Array<double> AP_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.dim();
  const ONE = 1.0, ZERO = 0.0;
  bool NOUNIT, UPPER;
  int J, JC, JCLAST = 0, JJ;
  double AJJ;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  NOUNIT = lsame(DIAG, 'N');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (!NOUNIT && !lsame(DIAG, 'U')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  }
  if (INFO.value != 0) {
    xerbla('DTPTRI', -INFO.value);
    return;
  }

  // Check for singularity if non-unit.

  if (NOUNIT) {
    if (UPPER) {
      JJ = 0;
      for (INFO.value = 1; INFO.value <= N; INFO.value++) {
        // 10
        JJ = JJ + INFO.value;
        if (AP[JJ] == ZERO) return;
      } // 10
    } else {
      JJ = 1;
      for (INFO.value = 1; INFO.value <= N; INFO.value++) {
        // 20
        if (AP[JJ] == ZERO) return;
        JJ = JJ + N - INFO.value + 1;
      } // 20
    }
    INFO.value = 0;
  }

  if (UPPER) {
    // Compute inverse of upper triangular matrix.

    JC = 1;
    for (J = 1; J <= N; J++) {
      // 30
      if (NOUNIT) {
        AP[JC + J - 1] = ONE / AP[JC + J - 1];
        AJJ = -AP[JC + J - 1];
      } else {
        AJJ = -ONE;
      }

      // Compute elements 1:j-1 of j-th column.

      dtpmv('Upper', 'No transpose', DIAG, J - 1, AP, AP(JC), 1);
      dscal(J - 1, AJJ, AP(JC), 1);
      JC = JC + J;
    } // 30
  } else {
    // Compute inverse of lower triangular matrix.

    JC = N * (N + 1) ~/ 2;
    for (J = N; J >= 1; J--) {
      // 40
      if (NOUNIT) {
        AP[JC] = ONE / AP[JC];
        AJJ = -AP[JC];
      } else {
        AJJ = -ONE;
      }
      if (J < N) {
        // Compute elements j+1:n of j-th column.

        dtpmv('Lower', 'No transpose', DIAG, N - J, AP(JCLAST), AP(JC + 1), 1);
        dscal(N - J, AJJ, AP(JC + 1), 1);
      }
      JCLAST = JC;
      JC = JC - N + J - 2;
    } // 40
  }
}
