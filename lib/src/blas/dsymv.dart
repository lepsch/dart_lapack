import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/matrix.dart';

void dsymv(
  final String UPLO,
  final int N,
  final double ALPHA,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> X_,
  final int INCX,
  final double BETA,
  final Array<double> Y_,
  final int INCY,
) {
// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final X = X_.dim();
  final Y = Y_.dim();
  const ONE = 1.0, ZERO = 0.0;
  double TEMP1, TEMP2;
  int I, INFO, IX, IY, J, JX, JY, KX, KY;

  // Test the input parameters.

  INFO = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO = 1;
  } else if (N < 0) {
    INFO = 2;
  } else if (LDA < max(1, N)) {
    INFO = 5;
  } else if (INCX == 0) {
    INFO = 7;
  } else if (INCY == 0) {
    INFO = 10;
  }
  if (INFO != 0) {
    xerbla('DSYMV ', INFO);
    return;
  }

  // Quick return if possible.

  if ((N == 0) || ((ALPHA == ZERO) && (BETA == ONE))) return;

  // Set up the start points in  X  and  Y.

  if (INCX > 0) {
    KX = 1;
  } else {
    KX = 1 - (N - 1) * INCX;
  }
  if (INCY > 0) {
    KY = 1;
  } else {
    KY = 1 - (N - 1) * INCY;
  }

  // Start the operations. In this version the elements of A are
  // accessed sequentially with one pass through the triangular part
  // of A.

  // First form  y := beta*y.

  if (BETA != ONE) {
    if (INCY == 1) {
      if (BETA == ZERO) {
        for (I = 1; I <= N; I++) {
          Y[I] = ZERO;
        }
      } else {
        for (I = 1; I <= N; I++) {
          Y[I] = BETA * Y[I];
        }
      }
    } else {
      IY = KY;
      if (BETA == ZERO) {
        for (I = 1; I <= N; I++) {
          Y[IY] = ZERO;
          IY = IY + INCY;
        }
      } else {
        for (I = 1; I <= N; I++) {
          Y[IY] = BETA * Y[IY];
          IY = IY + INCY;
        }
      }
    }
  }
  if (ALPHA == ZERO) return;
  if (lsame(UPLO, 'U')) {
    // Form  y  when A is stored in upper triangle.

    if ((INCX == 1) && (INCY == 1)) {
      for (J = 1; J <= N; J++) {
        TEMP1 = ALPHA * X[J];
        TEMP2 = ZERO;
        for (I = 1; I <= J - 1; I++) {
          Y[I] = Y[I] + TEMP1 * A[I][J];
          TEMP2 = TEMP2 + A[I][J] * X[I];
        }
        Y[J] = Y[J] + TEMP1 * A[J][J] + ALPHA * TEMP2;
      }
    } else {
      JX = KX;
      JY = KY;
      for (J = 1; J <= N; J++) {
        TEMP1 = ALPHA * X[JX];
        TEMP2 = ZERO;
        IX = KX;
        IY = KY;
        for (I = 1; I <= J - 1; I++) {
          Y[IY] = Y[IY] + TEMP1 * A[I][J];
          TEMP2 = TEMP2 + A[I][J] * X[IX];
          IX = IX + INCX;
          IY = IY + INCY;
        }
        Y[JY] = Y[JY] + TEMP1 * A[J][J] + ALPHA * TEMP2;
        JX = JX + INCX;
        JY = JY + INCY;
      }
    }
  } else {
    // Form  y  when A is stored in lower triangle.

    if ((INCX == 1) && (INCY == 1)) {
      for (J = 1; J <= N; J++) {
        TEMP1 = ALPHA * X[J];
        TEMP2 = ZERO;
        Y[J] = Y[J] + TEMP1 * A[J][J];
        for (I = J + 1; I <= N; I++) {
          Y[I] = Y[I] + TEMP1 * A[I][J];
          TEMP2 = TEMP2 + A[I][J] * X[I];
        }
        Y[J] = Y[J] + ALPHA * TEMP2;
      }
    } else {
      JX = KX;
      JY = KY;
      for (J = 1; J <= N; J++) {
        TEMP1 = ALPHA * X[JX];
        TEMP2 = ZERO;
        Y[JY] = Y[JY] + TEMP1 * A[J][J];
        IX = JX;
        IY = JY;
        for (I = J + 1; I <= N; I++) {
          IX = IX + INCX;
          IY = IY + INCY;
          Y[IY] = Y[IY] + TEMP1 * A[I][J];
          TEMP2 = TEMP2 + A[I][J] * X[IX];
        }
        Y[JY] = Y[JY] + ALPHA * TEMP2;
        JX = JX + INCX;
        JY = JY + INCY;
      }
    }
  }
}
