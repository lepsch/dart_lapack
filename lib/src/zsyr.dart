import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zsyr(
  final String UPLO,
  final int N,
  final Complex ALPHA,
  final Array<Complex> X_,
  final int INCX,
  final Matrix<Complex> A_,
  final int LDA,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having();
  final A = A_.having(ld: LDA);
  int I, INFO, IX, J, JX, KX = 0;
  Complex TEMP;

  // Test the input parameters.

  INFO = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO = 1;
  } else if (N < 0) {
    INFO = 2;
  } else if (INCX == 0) {
    INFO = 5;
  } else if (LDA < max(1, N)) {
    INFO = 7;
  }
  if (INFO != 0) {
    xerbla('ZSYR  ', INFO);
    return;
  }

  // Quick return if possible.

  if ((N == 0) || (ALPHA == Complex.zero)) return;

  // Set the start point in X if the increment is not unity.

  if (INCX <= 0) {
    KX = 1 - (N - 1) * INCX;
  } else if (INCX != 1) {
    KX = 1;
  }

  // Start the operations. In this version the elements of A are
  // accessed sequentially with one pass through the triangular part
  // of A.

  if (lsame(UPLO, 'U')) {
    // Form  A  when A is stored in upper triangle.

    if (INCX == 1) {
      for (J = 1; J <= N; J++) {
        // 20
        if (X[J] != Complex.zero) {
          TEMP = ALPHA * X[J];
          for (I = 1; I <= J; I++) {
            // 10
            A[I][J] = A[I][J] + X[I] * TEMP;
          } // 10
        }
      } // 20
    } else {
      JX = KX;
      for (J = 1; J <= N; J++) {
        // 40
        if (X[JX] != Complex.zero) {
          TEMP = ALPHA * X[JX];
          IX = KX;
          for (I = 1; I <= J; I++) {
            // 30
            A[I][J] = A[I][J] + X[IX] * TEMP;
            IX = IX + INCX;
          } // 30
        }
        JX = JX + INCX;
      } // 40
    }
  } else {
    // Form  A  when A is stored in lower triangle.

    if (INCX == 1) {
      for (J = 1; J <= N; J++) {
        // 60
        if (X[J] != Complex.zero) {
          TEMP = ALPHA * X[J];
          for (I = J; I <= N; I++) {
            // 50
            A[I][J] = A[I][J] + X[I] * TEMP;
          } // 50
        }
      } // 60
    } else {
      JX = KX;
      for (J = 1; J <= N; J++) {
        // 80
        if (X[JX] != Complex.zero) {
          TEMP = ALPHA * X[JX];
          IX = JX;
          for (I = J; I <= N; I++) {
            // 70
            A[I][J] = A[I][J] + X[IX] * TEMP;
            IX = IX + INCX;
          } // 70
        }
        JX = JX + INCX;
      } // 80
    }
  }
}
