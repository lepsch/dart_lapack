import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/matrix.dart';

void dsyr2(
  final String UPLO,
  final int N,
  final double ALPHA,
  final Array<double> X_,
  final int INCX,
  final Array<double> Y_,
  final int INCY,
  final Matrix<double> A_,
  final int LDA,
) {
// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having();
  final Y = Y_.having();
  final A = A_.having(ld: LDA);
  const ZERO = 0.0;

  // Test the input parameters.

  var INFO = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO = 1;
  } else if (N < 0) {
    INFO = 2;
  } else if (INCX == 0) {
    INFO = 5;
  } else if (INCY == 0) {
    INFO = 7;
  } else if (LDA < max(1, N)) {
    INFO = 9;
  }
  if (INFO != 0) {
    xerbla('DSYR2 ', INFO);
    return;
  }

  // Quick return if possible.

  if ((N == 0) || (ALPHA == ZERO)) return;

  // Set up the start points in X and Y if the increments are not both
  // unity.

  int JX = 0, JY = 0, KX = 0, KY = 0;
  if ((INCX != 1) || (INCY != 1)) {
    KX = INCX > 0 ? 1 : 1 - (N - 1) * INCX;
    KY = INCY > 0 ? 1 : 1 - (N - 1) * INCY;
    JX = KX;
    JY = KY;
  }

  // Start the operations. In this version the elements of A are
  // accessed sequentially with one pass through the triangular part
  // of A.

  if (lsame(UPLO, 'U')) {
    // Form  A  when A is stored in the upper triangle.

    if ((INCX == 1) && (INCY == 1)) {
      for (var J = 1; J <= N; J++) {
        if ((X[J] != ZERO) || (Y[J] != ZERO)) {
          final TEMP1 = ALPHA * Y[J];
          final TEMP2 = ALPHA * X[J];
          for (var I = 1; I <= J; I++) {
            A[I][J] += X[I] * TEMP1 + Y[I] * TEMP2;
          }
        }
      }
    } else {
      for (var J = 1; J <= N; J++) {
        if ((X[JX] != ZERO) || (Y[JY] != ZERO)) {
          final TEMP1 = ALPHA * Y[JY];
          final TEMP2 = ALPHA * X[JX];
          var IX = KX, IY = KY;
          for (var I = 1; I <= J; I++) {
            A[I][J] += X[IX] * TEMP1 + Y[IY] * TEMP2;
            IX += INCX;
            IY += INCY;
          }
        }
        JX += INCX;
        JY += INCY;
      }
    }
  } else {
    // Form  A  when A is stored in the lower triangle.

    if ((INCX == 1) && (INCY == 1)) {
      for (var J = 1; J <= N; J++) {
        if ((X[J] != ZERO) || (Y[J] != ZERO)) {
          final TEMP1 = ALPHA * Y[J];
          final TEMP2 = ALPHA * X[J];
          for (var I = J; I <= N; I++) {
            A[I][J] += X[I] * TEMP1 + Y[I] * TEMP2;
          }
        }
      }
    } else {
      for (var J = 1; J <= N; J++) {
        if ((X[JX] != ZERO) || (Y[JY] != ZERO)) {
          final TEMP1 = ALPHA * Y[JY];
          final TEMP2 = ALPHA * X[JX];
          var IX = JX, IY = JY;
          for (var I = J; I <= N; I++) {
            A[I][J] += X[IX] * TEMP1 + Y[IY] * TEMP2;
            IX += INCX;
            IY += INCY;
          }
        }
        JX += INCX;
        JY += INCY;
      }
    }
  }
}
