import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/matrix.dart';

void dspr2(
  final String UPLO,
  final int N,
  final double ALPHA,
  final Array<double> X,
  final int INCX,
  final Array<double> Y,
  final int INCY,
  final Array<double> AP,
) {
// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;
  double TEMP1, TEMP2;
  int I, INFO, IX, IY, J, JX = 0, JY = 0, K, KK, KX = 0, KY = 0;

  // Test the input parameters.

  INFO = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO = 1;
  } else if (N < 0) {
    INFO = 2;
  } else if (INCX == 0) {
    INFO = 5;
  } else if (INCY == 0) {
    INFO = 7;
  }
  if (INFO != 0) {
    xerbla('DSPR2 ', INFO);
    return;
  }

  // Quick return if possible.

  if ((N == 0) || (ALPHA == ZERO)) return;

  // Set up the start points in X and Y if the increments are not both
  // unity.

  if ((INCX != 1) || (INCY != 1)) {
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
    JX = KX;
    JY = KY;
  }

  // Start the operations. In this version the elements of the array AP
  // are accessed sequentially with one pass through AP.

  KK = 1;
  if (lsame(UPLO, 'U')) {
    // Form  A  when upper triangle is stored in AP.

    if ((INCX == 1) && (INCY == 1)) {
      for (J = 1; J <= N; J++) {
        if ((X[J] != ZERO) || (Y[J] != ZERO)) {
          TEMP1 = ALPHA * Y[J];
          TEMP2 = ALPHA * X[J];
          K = KK;
          for (I = 1; I <= J; I++) {
            AP[K] = AP[K] + X[I] * TEMP1 + Y[I] * TEMP2;
            K = K + 1;
          }
        }
        KK = KK + J;
      }
    } else {
      for (J = 1; J <= N; J++) {
        if ((X[JX] != ZERO) || (Y[JY] != ZERO)) {
          TEMP1 = ALPHA * Y[JY];
          TEMP2 = ALPHA * X[JX];
          IX = KX;
          IY = KY;
          for (K = KK; K <= KK + J - 1; K++) {
            AP[K] = AP[K] + X[IX] * TEMP1 + Y[IY] * TEMP2;
            IX = IX + INCX;
            IY = IY + INCY;
          }
        }
        JX = JX + INCX;
        JY = JY + INCY;
        KK = KK + J;
      }
    }
  } else {
    // Form  A  when lower triangle is stored in AP.

    if ((INCX == 1) && (INCY == 1)) {
      for (J = 1; J <= N; J++) {
        if ((X[J] != ZERO) || (Y[J] != ZERO)) {
          TEMP1 = ALPHA * Y[J];
          TEMP2 = ALPHA * X[J];
          K = KK;
          for (I = J; I <= N; I++) {
            AP[K] = AP[K] + X[I] * TEMP1 + Y[I] * TEMP2;
            K = K + 1;
          }
        }
        KK = KK + N - J + 1;
      }
    } else {
      for (J = 1; J <= N; J++) {
        if ((X[JX] != ZERO) || (Y[JY] != ZERO)) {
          TEMP1 = ALPHA * Y[JY];
          TEMP2 = ALPHA * X[JX];
          IX = JX;
          IY = JY;
          for (K = KK; K <= KK + N - J; K++) {
            AP[K] = AP[K] + X[IX] * TEMP1 + Y[IY] * TEMP2;
            IX = IX + INCX;
            IY = IY + INCY;
          }
        }
        JX = JX + INCX;
        JY = JY + INCY;
        KK = KK + N - J + 1;
      }
    }
  }
}
