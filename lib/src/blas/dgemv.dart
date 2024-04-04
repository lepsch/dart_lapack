import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/matrix.dart';

void dgemv(
  final String TRANS,
  final int M,
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
  final A = A_.having(ld: LDA);
  final X = X_.having();
  final Y = Y_.having();
  const ONE = 1.0, ZERO = 0.0;

  // Test the input parameters.

  var INFO = 0;
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
    xerbla('DGEMV', INFO);
    return;
  }

  // Quick return if possible.

  if ((M == 0) || (N == 0) || ((ALPHA == ZERO) && (BETA == ONE))) return;

  // Set  LENX  and  LENY, the lengths of the vectors x and y, and set
  // up the start points in  X  and  Y.

  final (LENX, LENY) = lsame(TRANS, 'N') ? (N, M) : (M, N);
  final KX = INCX > 0 ? 1 : 1 - (LENX - 1) * INCX;
  final KY = INCY > 0 ? 1 : 1 - (LENY - 1) * INCY;

  // Start the operations. In this version the elements of A are
  // accessed sequentially with one pass through A.

  // First form  y := beta*y.

  if (BETA != ONE) {
    if (INCY == 1) {
      if (BETA == ZERO) {
        for (var I = 1; I <= LENY; I++) {
          Y[I] = ZERO;
        }
      } else {
        for (var I = 1; I <= LENY; I++) {
          Y[I] *= BETA;
        }
      }
    } else {
      var IY = KY;
      if (BETA == ZERO) {
        for (var I = 1; I <= LENY; I++) {
          Y[IY] = ZERO;
          IY += INCY;
        }
      } else {
        for (var I = 1; I <= LENY; I++) {
          Y[IY] *= BETA;
          IY += INCY;
        }
      }
    }
  }
  if (ALPHA == ZERO) return;
  if (lsame(TRANS, 'N')) {
    // Form  y := alpha*A*x + y.

    var JX = KX;
    if (INCY == 1) {
      for (var J = 1; J <= N; J++) {
        final TEMP = ALPHA * X[JX];
        for (var I = 1; I <= M; I++) {
          Y[I] += TEMP * A[I][J];
        }
        JX += INCX;
      }
    } else {
      for (var J = 1; J <= N; J++) {
        final TEMP = ALPHA * X[JX];
        var IY = KY;
        for (var I = 1; I <= M; I++) {
          Y[IY] += TEMP * A[I][J];
          IY += INCY;
        }
        JX += INCX;
      }
    }
  } else {
    // Form  y := alpha*A**T*x + y.
    var JY = KY;
    if (INCX == 1) {
      for (var J = 1; J <= N; J++) {
        var TEMP = ZERO;
        for (var I = 1; I <= M; I++) {
          TEMP += A[I][J] * X[I];
        }
        Y[JY] += ALPHA * TEMP;
        JY += INCY;
      }
    } else {
      for (var J = 1; J <= N; J++) {
        var TEMP = ZERO;
        var IX = KX;
        for (var I = 1; I <= M; I++) {
          TEMP += A[I][J] * X[IX];
          IX += INCX;
        }
        Y[JY] += ALPHA * TEMP;
        JY += INCY;
      }
    }
  }
}
