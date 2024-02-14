import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zgbmv(
  final String TRANS,
  final int M,
  final int N,
  final int KL,
  final int KU,
  final double ALPHA,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> X_,
  final int INCX,
  final double BETA,
  final Array<Complex> Y_,
  final int INCY,
) {
// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final X = X_.dim();
  final Y = Y_.dim();
  Complex TEMP;
  int I, INFO, IX, IY, J, JX, JY, K, KUP1, KX, KY, LENX, LENY;
  bool NOCONJ;

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
    xerbla('ZGBMV ', INFO);
    return;
  }

  // Quick return if possible.

  if ((M == 0) ||
      (N == 0) ||
      ((ALPHA.toComplex() == Complex.zero) &&
          (BETA.toComplex() == Complex.one))) return;

  NOCONJ = lsame(TRANS, 'T');

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

  if (BETA.toComplex() != Complex.one) {
    if (INCY == 1) {
      if (BETA.toComplex() == Complex.zero) {
        for (I = 1; I <= LENY; I++) {
          Y[I] = Complex.zero;
        }
      } else {
        for (I = 1; I <= LENY; I++) {
          Y[I] = Complex(BETA) * Y[I];
        }
      }
    } else {
      IY = KY;
      if (BETA.toComplex() == Complex.zero) {
        for (I = 1; I <= LENY; I++) {
          Y[IY] = Complex.zero;
          IY = IY + INCY;
        }
      } else {
        for (I = 1; I <= LENY; I++) {
          Y[IY] = Complex(BETA) * Y[IY];
          IY = IY + INCY;
        }
      }
    }
  }
  if (ALPHA.toComplex() == Complex.zero) return;
  KUP1 = KU + 1;
  if (lsame(TRANS, 'N')) {
    // Form  y := alpha*A*x + y.

    JX = KX;
    if (INCY == 1) {
      for (J = 1; J <= N; J++) {
        TEMP = Complex(ALPHA) * X[JX];
        K = KUP1 - J;
        for (I = max(1, J - KU); I <= min(M, J + KL); I++) {
          Y[I] = Y[I] + TEMP * A[K + I][J];
        }
        JX = JX + INCX;
      }
    } else {
      for (J = 1; J <= N; J++) {
        TEMP = Complex(ALPHA) * X[JX];
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
    // Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y.

    JY = KY;
    if (INCX == 1) {
      for (J = 1; J <= N; J++) {
        TEMP = Complex.zero;
        K = KUP1 - J;
        if (NOCONJ) {
          for (I = max(1, J - KU); I <= min(M, J + KL); I++) {
            TEMP = TEMP + A[K + I][J] * X[I];
          }
        } else {
          for (I = max(1, J - KU); I <= min(M, J + KL); I++) {
            TEMP = TEMP + A[K + I][J].conjugate() * X[I];
          }
        }
        Y[JY] = Y[JY] + Complex(ALPHA) * TEMP;
        JY = JY + INCY;
      }
    } else {
      for (J = 1; J <= N; J++) {
        TEMP = Complex.zero;
        IX = KX;
        K = KUP1 - J;
        if (NOCONJ) {
          for (I = max(1, J - KU); I <= min(M, J + KL); I++) {
            TEMP = TEMP + A[K + I][J] * X[IX];
            IX = IX + INCX;
          }
        } else {
          for (I = max(1, J - KU); I <= min(M, J + KL); I++) {
            TEMP = TEMP + (A[K + I][J]).conjugate() * X[IX];
            IX = IX + INCX;
          }
        }
        Y[JY] = Y[JY] + Complex(ALPHA) * TEMP;
        JY = JY + INCY;
        if (J > KU) KX = KX + INCX;
      }
    }
  }
}
