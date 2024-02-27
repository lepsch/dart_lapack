import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zdotc.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacgv.dart';

void zpotf2(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  const ONE = 1.0, ZERO = 0.0;
  bool UPPER;
  int J;
  double AJJ;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('ZPOTF2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (UPPER) {
    // Compute the Cholesky factorization A = U**H *U.

    for (J = 1; J <= N; J++) {
      // 10

      // Compute U(J,J) and test for non-positive-definiteness.

      AJJ = A[J][J].toDouble() -
          (zdotc(J - 1, A(1, J).asArray(), 1, A(1, J).asArray(), 1)).toDouble();
      if (AJJ <= ZERO || disnan(AJJ)) {
        A[J][J] = AJJ.toComplex();
        INFO.value = J;
        return;
      }
      AJJ = sqrt(AJJ);
      A[J][J] = AJJ.toComplex();

      // Compute elements J+1:N of row J.

      if (J < N) {
        zlacgv(J - 1, A(1, J).asArray(), 1);
        zgemv('Transpose', J - 1, N - J, -Complex.one, A(1, J + 1), LDA,
            A(1, J).asArray(), 1, Complex.one, A(J, J + 1).asArray(), LDA);
        zlacgv(J - 1, A(1, J).asArray(), 1);
        zdscal(N - J, ONE / AJJ, A(J, J + 1).asArray(), LDA);
      }
    } // 10
  } else {
    // Compute the Cholesky factorization A = L*L**H.

    for (J = 1; J <= N; J++) {
      // 20

      // Compute L(J,J) and test for non-positive-definiteness.

      AJJ = A[J][J].toDouble() -
          (zdotc(J - 1, A(J, 1).asArray(), LDA, A(J, 1).asArray(), LDA))
              .toDouble();
      if (AJJ <= ZERO || disnan(AJJ)) {
        A[J][J] = AJJ.toComplex();
        INFO.value = J;
        return;
      }
      AJJ = sqrt(AJJ);
      A[J][J] = AJJ.toComplex();

      // Compute elements J+1:N of column J.

      if (J < N) {
        zlacgv(J - 1, A(J, 1).asArray(), LDA);
        zgemv('No transpose', N - J, J - 1, -Complex.one, A(J + 1, 1), LDA,
            A(J, 1).asArray(), LDA, Complex.one, A(J + 1, J).asArray(), 1);
        zlacgv(J - 1, A(J, 1).asArray(), LDA);
        zdscal(N - J, ONE / AJJ, A(J + 1, J).asArray(), 1);
      }
    } // 20
  }
}
