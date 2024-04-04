import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zgemv(
  final String TRANS,
  final int M,
  final int N,
  final Complex ALPHA,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> X_,
  final int INCX,
  final Complex BETA,
  final Array<Complex> Y_,
  final int INCY,
) {
// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final X = X_.having();
  final Y = Y_.having();
  Complex TEMP;
  int I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY;
  bool NOCONJ;

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
    xerbla('ZGEMV', INFO);
    return;
  }

  // Quick return if possible.

  if ((M == 0) ||
      (N == 0) ||
      ((ALPHA == Complex.zero) && (BETA == Complex.one))) return;

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
  // accessed sequentially with one pass through A.

  // First form  y := beta*y.

  if (BETA != Complex.one) {
    if (INCY == 1) {
      if (BETA == Complex.zero) {
        for (I = 1; I <= LENY; I++) {
          Y[I] = Complex.zero;
        }
      } else {
        for (I = 1; I <= LENY; I++) {
          Y[I] = BETA * Y[I];
        }
      }
    } else {
      IY = KY;
      if (BETA == Complex.zero) {
        for (I = 1; I <= LENY; I++) {
          Y[IY] = Complex.zero;
          IY += INCY;
        }
      } else {
        for (I = 1; I <= LENY; I++) {
          Y[IY] = BETA * Y[IY];
          IY += INCY;
        }
      }
    }
  }
  if (ALPHA == Complex.zero) return;
  if (lsame(TRANS, 'N')) {
    // Form  y := alpha*A*x + y.

    JX = KX;
    if (INCY == 1) {
      for (J = 1; J <= N; J++) {
        TEMP = ALPHA * X[JX];
        for (I = 1; I <= M; I++) {
          Y[I] += TEMP * A[I][J];
        }
        JX += INCX;
      }
    } else {
      for (J = 1; J <= N; J++) {
        TEMP = ALPHA * X[JX];
        IY = KY;
        for (I = 1; I <= M; I++) {
          Y[IY] += TEMP * A[I][J];
          IY += INCY;
        }
        JX += INCX;
      }
    }
  } else {
    // Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y.

    JY = KY;
    if (INCX == 1) {
      for (J = 1; J <= N; J++) {
        TEMP = Complex.zero;
        if (NOCONJ) {
          for (I = 1; I <= M; I++) {
            TEMP += A[I][J] * X[I];
          }
        } else {
          for (I = 1; I <= M; I++) {
            TEMP += A[I][J].conjugate() * X[I];
          }
        }
        Y[JY] += ALPHA * TEMP;
        JY += INCY;
      }
    } else {
      for (J = 1; J <= N; J++) {
        TEMP = Complex.zero;
        IX = KX;
        if (NOCONJ) {
          for (I = 1; I <= M; I++) {
            TEMP += A[I][J] * X[IX];
            IX += INCX;
          }
        } else {
          for (I = 1; I <= M; I++) {
            TEMP += A[I][J].conjugate() * X[IX];
            IX += INCX;
          }
        }
        Y[JY] += ALPHA * TEMP;
        JY += INCY;
      }
    }
  }
}
