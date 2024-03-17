import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dspr.dart';
import 'package:lapack/src/blas/dtpmv.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dtptri.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dpptri(
  final String UPLO,
  final int N,
  final Array<double> AP_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.having();
  const ONE = 1.0;
  bool UPPER;
  int J, JC, JJ, JJN;
  double AJJ;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  }
  if (INFO.value != 0) {
    xerbla('DPPTRI', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Invert the triangular Cholesky factor U or L.

  dtptri(UPLO, 'Non-unit', N, AP, INFO);
  if (INFO.value > 0) return;

  if (UPPER) {
    // Compute the product inv(U) * inv(U)**T.

    JJ = 0;
    for (J = 1; J <= N; J++) {
      JC = JJ + 1;
      JJ += J;
      if (J > 1) dspr('Upper', J - 1, ONE, AP(JC), 1, AP);
      AJJ = AP[JJ];
      dscal(J, AJJ, AP(JC), 1);
    }
  } else {
    // Compute the product inv(L)**T * inv(L).

    JJ = 1;
    for (J = 1; J <= N; J++) {
      JJN = JJ + N - J + 1;
      AP[JJ] = ddot(N - J + 1, AP(JJ), 1, AP(JJ), 1);
      if (J < N) {
        dtpmv('Lower', 'Transpose', 'Non-unit', N - J, AP(JJN), AP(JJ + 1), 1);
      }
      JJ = JJN;
    }
  }
}
