import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/matrix.dart';

void dgbmv(
  final String TRANS,
  final int M,
  final int N,
  final int KL,
  final int KU,
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
  double TEMP;
  int I, INFO, IX, IY, J, JX, JY, K, KUP1, KX, KY, LENX, LENY;

  // Test the input parameters.

  INFO = 0;
  if (!lsame(TRANS, 'N') && !lsame(TRANS, 'T') && !lsame(TRANS, 'C')) {
    INFO = 1;
  } else if (M < 0) {
    INFO = 2;
  } else if (N < 0) {
    INFO = 3;
  } else if (KL < 0) {
    INFO = 4;
  } else if (KU < 0) {
    INFO = 5;
  } else if (LDA < (KL + KU + 1)) {
    INFO = 8;
  } else if (INCX == 0) {
    INFO = 10;
  } else if (INCY == 0) {
    INFO = 13;
  }
  if (INFO != 0) {
    xerbla('DGBMV ', INFO);
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
  // accessed sequentially with one pass through the band part of A.

  // First form  y := beta*y.

  if (BETA != ONE) {
    if (INCY == 1) {
      if (BETA == ZERO) {
        for (I = 1; I <= LENY; I++) {
          Y[I] = ZERO;
        }
      } else {
        for (I = 1; I <= LENY; I++) {
          Y[I] = BETA * Y[I];
        }
      }
    } else {
      IY = KY;
      if (BETA == ZERO) {
        for (I = 1; I <= LENY; I++) {
          Y[IY] = ZERO;
          IY = IY + INCY;
        }
      } else {
        for (I = 1; I <= LENY; I++) {
          Y[IY] = BETA * Y[IY];
          IY = IY + INCY;
        }
      }
    }
  }
  if (ALPHA == ZERO) return;
  KUP1 = KU + 1;
  if (lsame(TRANS, 'N')) {
    // Form  y := alpha*A*x + y.
    JX = KX;
    if (INCY == 1) {
      for (J = 1; J <= N; J++) {
        TEMP = ALPHA * X[JX];
        K = KUP1 - J;
        for (I = max(1, J - KU); I <= min(M, J + KL); I++) {
          Y[I] = Y[I] + TEMP * A[K + I][J];
        }
        JX = JX + INCX;
      }
    } else {
      for (J = 1; J <= N; J++) {
        TEMP = ALPHA * X[JX];
        IY = KY;
        K = KUP1 - J;
        for (I = max(1, J - KU); I <= min(M, J + KL); I++) {
          Y[IY] = Y[IY] + TEMP * A[K + I][J];
          IY = IY + INCY;
        }
        JX = JX + INCX;
        if (J > KU) KY = KY + INCY;
      }
    }
  } else {
    // Form  y := alpha*A**T*x + y.

    JY = KY;
    if (INCX == 1) {
      for (J = 1; J <= N; J++) {
        TEMP = ZERO;
        K = KUP1 - J;
        for (I = max(1, J - KU); I <= min(M, J + KL); I++) {
          TEMP = TEMP + A[K + I][J] * X[I];
        }
        Y[JY] = Y[JY] + ALPHA * TEMP;
        JY = JY + INCY;
      }
    } else {
      for (J = 1; J <= N; J++) {
        TEMP = ZERO;
        IX = KX;
        K = KUP1 - J;
        for (I = max(1, J - KU); I <= min(M, J + KL); I++) {
          TEMP = TEMP + A[K + I][J] * X[IX];
          IX = IX + INCX;
        }
        Y[JY] = Y[JY] + ALPHA * TEMP;
        JY = JY + INCY;
        if (J > KU) KX = KX + INCX;
      }
    }
  }
}
