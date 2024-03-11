import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zdotc.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zhpr.dart';
import 'package:lapack/src/blas/ztpmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/ztptri.dart';

void zpptri(
  final String UPLO,
  final int N,
  final Array<Complex> AP_,
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
    xerbla('ZPPTRI', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Invert the triangular Cholesky factor U or L.

  ztptri(UPLO, 'Non-unit', N, AP, INFO);
  if (INFO.value > 0) return;
  if (UPPER) {
    // Compute the product inv(U) * inv(U)**H.

    JJ = 0;
    for (J = 1; J <= N; J++) {
      // 10
      JC = JJ + 1;
      JJ += J;
      if (J > 1) zhpr('Upper', J - 1, ONE, AP(JC), 1, AP);
      AJJ = (AP[JJ]).toDouble();
      zdscal(J, AJJ, AP(JC), 1);
    } // 10
  } else {
    // Compute the product inv(L)**H * inv(L).

    JJ = 1;
    for (J = 1; J <= N; J++) {
      // 20
      JJN = JJ + N - J + 1;
      AP[JJ] = zdotc(N - J + 1, AP(JJ), 1, AP(JJ), 1).real.toComplex();
      if (J < N) {
        ztpmv('Lower', 'Conjugate transpose', 'Non-unit', N - J, AP(JJN),
            AP(JJ + 1), 1);
      }
      JJ = JJN;
    } // 20
  }
}
