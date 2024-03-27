import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zdotc.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zhpr.dart';
import 'package:lapack/src/blas/ztpsv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zpptrf(
  final String UPLO,
  final int N,
  final Array<Complex> AP_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool UPPER;
  int J, JC, JJ;
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
    xerbla('ZPPTRF', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (UPPER) {
    // Compute the Cholesky factorization A = U**H * U.

    JJ = 0;
    for (J = 1; J <= N; J++) {
      JC = JJ + 1;
      JJ += J;

      // Compute elements 1:J-1 of column J.

      if (J > 1) {
        ztpsv('Upper', 'Conjugate transpose', 'Non-unit', J - 1, AP, AP(JC), 1);
      }

      // Compute U(J,J) and test for non-positive-definiteness.

      AJJ = AP[JJ].real - zdotc(J - 1, AP(JC), 1, AP(JC), 1).real;
      if (AJJ <= ZERO) {
        AP[JJ] = AJJ.toComplex();
        INFO.value = J;
        return;
      }
      AP[JJ] = sqrt(AJJ).toComplex();
    }
  } else {
    // Compute the Cholesky factorization A = L * L**H.

    JJ = 1;
    for (J = 1; J <= N; J++) {
      // Compute L(J,J) and test for non-positive-definiteness.

      AJJ = AP[JJ].real;
      if (AJJ <= ZERO) {
        AP[JJ] = AJJ.toComplex();
        INFO.value = J;
        return;
      }
      AJJ = sqrt(AJJ);
      AP[JJ] = AJJ.toComplex();

      // Compute elements J+1:N of column J and update the trailing
      // submatrix.

      if (J < N) {
        zdscal(N - J, ONE / AJJ, AP(JJ + 1), 1);
        zhpr('Lower', N - J, -ONE, AP(JJ + 1), 1, AP(JJ + N - J + 1));
        JJ += N - J + 1;
      }
    }
  }
}
