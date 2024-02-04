import 'dart:math';

import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/matrix.dart';

void dgemv(
  final String TRANS,
  final int M,
  final int N,
  final double ALPHA,
  final Matrix2d<double> A,
  final int LDA,
  final Array<double> X,
  final int INCX,
  final double BETA,
  final Array<double> Y,
  final int INCY,
) {
// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0, ZERO = 0.0;
  double TEMP;
  int I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY;

  // Test the input parameters.

  INFO = 0;
  if (!lsame(TRANS, 'N') && !lsame(TRANS, 'T') && !lsame(TRANS, 'C')) {
    INFO = 1;
  } else if (M < 0) {
    INFO = 2;
  } else if (N < 0) {
    INFO = 3;
  } else if (LDA < max(1, M)) {
    INFO = 6;
  } else if (INCX == 0) {
    INFO = 8;
  } else if (INCY == 0) {
    INFO = 11;
  }
  if (INFO != 0) {
    xerbla('DGEMV ', INFO);
    return;
  }

  // Quick return if possible.

  if ((M == 0) || (N == 0) || ((ALPHA == ZERO) && (BETA == ONE))) return;

  // Set  LENX  and  LENY, the lengths of the vectors x and y, and set
  // up the start points in  X  and  Y.

  if (lsame(TRANS, 'N')) {
    LENX = N;
    LENY = M;
  } else {
    LENX = M;
    LENY = N;
  }
  if (INCX > 0) {
    KX = 1;
  } else {
    KX = 1 - (LENX - 1) * INCX;
  }
  if (INCY > 0) {
    KY = 1;
  } else {
    KY = 1 - (LENY - 1) * INCY;
  }

  // Start the operations. In this version the elements of A are
  // accessed sequentially with one pass through A.

  // First form  y := beta*y.

  if (BETA != ONE) {
    if (INCY == 1) {
      if (BETA == ZERO) {
        for (I = 1; I <= LENY; I++) {
          // 10
          Y[I] = ZERO;
        } // 10
      } else {
        for (I = 1; I <= LENY; I++) {
          // 20
          Y[I] = BETA * Y(I);
        } // 20
      }
    } else {
      IY = KY;
      if (BETA == ZERO) {
        for (I = 1; I <= LENY; I++) {
          // 30
          Y[IY] = ZERO;
          IY = IY + INCY;
        } // 30
      } else {
        for (I = 1; I <= LENY; I++) {
          // 40
          Y[IY] = BETA * Y(IY);
          IY = IY + INCY;
        } // 40
      }
    }
  }
  if (ALPHA == ZERO) return;
  if (lsame(TRANS, 'N')) {
    // Form  y := alpha*A*x + y.

    JX = KX;
    if (INCY == 1) {
      for (J = 1; J <= N; J++) {
        // 60
        TEMP = ALPHA * X(JX);
        for (I = 1; I <= M; I++) {
          // 50
          Y[I] = Y(I) + TEMP * A(I, J);
        } // 50
        JX = JX + INCX;
      } // 60
    } else {
      for (J = 1; J <= N; J++) {
        // 80
        TEMP = ALPHA * X(JX);
        IY = KY;
        for (I = 1; I <= M; I++) {
          // 70
          Y[IY] = Y(IY) + TEMP * A(I, J);
          IY = IY + INCY;
        } // 70
        JX = JX + INCX;
      } // 80
    }
  } else {
    // Form  y := alpha*A**T*x + y.

    JY = KY;
    if (INCX == 1) {
      for (J = 1; J <= N; J++) {
        // 100
        TEMP = ZERO;
        for (I = 1; I <= M; I++) {
          // 90
          TEMP = TEMP + A(I, J) * X(I);
        } // 90
        Y[JY] = Y(JY) + ALPHA * TEMP;
        JY = JY + INCY;
      } // 100
    } else {
      for (J = 1; J <= N; J++) {
        // 120
        TEMP = ZERO;
        IX = KX;
        for (I = 1; I <= M; I++) {
          // 110
          TEMP = TEMP + A(I, J) * X(IX);
          IX = IX + INCX;
        } // 110
        Y[JY] = Y(JY) + ALPHA * TEMP;
        JY = JY + INCY;
      } // 120
    }
  }

  return;
}
