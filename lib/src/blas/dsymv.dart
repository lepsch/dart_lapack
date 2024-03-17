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
  final A = A_.having(ld: LDA);
  final X = X_.having();
  final Y = Y_.having();
  const ONE = 1.0, ZERO = 0.0;

  // Test the input parameters.

  var INFO = 0;
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

  final KX = INCX > 0 ? 1 : 1 - (N - 1) * INCX;
  final KY = INCY > 0 ? 1 : 1 - (N - 1) * INCY;

  // Start the operations. In this version the elements of A are
  // accessed sequentially with one pass through the triangular part
  // of A.

  // First form  y := beta*y.

  if (BETA != ONE) {
    if (INCY == 1) {
      if (BETA == ZERO) {
        for (var I = 1; I <= N; I++) {
          Y[I] = ZERO;
        }
      } else {
        for (var I = 1; I <= N; I++) {
          Y[I] *= BETA;
        }
      }
    } else {
      var IY = KY;
      if (BETA == ZERO) {
        for (var I = 1; I <= N; I++) {
          Y[IY] = ZERO;
          IY += INCY;
        }
      } else {
        for (var I = 1; I <= N; I++) {
          Y[IY] *= BETA;
          IY += INCY;
        }
      }
    }
  }
  if (ALPHA == ZERO) return;
  if (lsame(UPLO, 'U')) {
    // Form  y  when A is stored in upper triangle.

    if ((INCX == 1) && (INCY == 1)) {
      for (var J = 1; J <= N; J++) {
        final TEMP1 = ALPHA * X[J];
        var TEMP2 = ZERO;
        for (var I = 1; I <= J - 1; I++) {
          Y[I] += TEMP1 * A[I][J];
          TEMP2 += A[I][J] * X[I];
        }
        Y[J] += TEMP1 * A[J][J] + ALPHA * TEMP2;
      }
    } else {
      var JX = KX;
      var JY = KY;
      for (var J = 1; J <= N; J++) {
        final TEMP1 = ALPHA * X[JX];
        var TEMP2 = ZERO;
        var IX = KX;
        var IY = KY;
        for (var I = 1; I <= J - 1; I++) {
          Y[IY] += TEMP1 * A[I][J];
          TEMP2 += A[I][J] * X[IX];
          IX += INCX;
          IY += INCY;
        }
        Y[JY] += TEMP1 * A[J][J] + ALPHA * TEMP2;
        JX += INCX;
        JY += INCY;
      }
    }
  } else {
    // Form  y  when A is stored in lower triangle.

    if ((INCX == 1) && (INCY == 1)) {
      for (var J = 1; J <= N; J++) {
        final TEMP1 = ALPHA * X[J];
        var TEMP2 = ZERO;
        Y[J] += TEMP1 * A[J][J];
        for (var I = J + 1; I <= N; I++) {
          Y[I] += TEMP1 * A[I][J];
          TEMP2 += A[I][J] * X[I];
        }
        Y[J] += ALPHA * TEMP2;
      }
    } else {
      var JX = KX;
      var JY = KY;
      for (var J = 1; J <= N; J++) {
        final TEMP1 = ALPHA * X[JX];
        var TEMP2 = ZERO;
        Y[JY] += TEMP1 * A[J][J];
        var IX = JX, IY = JY;
        for (var I = J + 1; I <= N; I++) {
          IX += INCX;
          IY += INCY;
          Y[IY] += TEMP1 * A[I][J];
          TEMP2 += A[I][J] * X[IX];
        }
        Y[JY] += ALPHA * TEMP2;
        JX += INCX;
        JY += INCY;
      }
    }
  }
}
