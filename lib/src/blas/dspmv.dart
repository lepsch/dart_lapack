import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/matrix.dart';

void dspmv(
  final String UPLO,
  final int N,
  final double ALPHA,
  final Array<double> AP_,
  final Array<double> X_,
  final int INCX,
  final double BETA,
  final Array<double> Y_,
  final int INCY,
) {
// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.dim();
  final X = X_.dim();
  final Y = Y_.dim();
  const ONE = 1.0, ZERO = 0.0;
  double TEMP1, TEMP2;
  int I, INFO, IX, IY, J, JX, JY, K, KK, KX, KY;

  // Test the input parameters.

  INFO = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO = 1;
  } else if (N < 0) {
    INFO = 2;
  } else if (INCX == 0) {
    INFO = 6;
  } else if (INCY == 0) {
    INFO = 9;
  }
  if (INFO != 0) {
    xerbla('DSPMV ', INFO);
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

  // Start the operations. In this version the elements of the array AP
  // are accessed sequentially with one pass through AP.

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
  KK = 1;
  if (lsame(UPLO, 'U')) {
    // Form  y  when AP contains the upper triangle.

    if ((INCX == 1) && (INCY == 1)) {
      for (J = 1; J <= N; J++) {
        TEMP1 = ALPHA * X[J];
        TEMP2 = ZERO;
        K = KK;
        for (I = 1; I <= J - 1; I++) {
          Y[I] = Y[I] + TEMP1 * AP[K];
          TEMP2 = TEMP2 + AP[K] * X[I];
          K = K + 1;
        }
        Y[J] = Y[J] + TEMP1 * AP[KK + J - 1] + ALPHA * TEMP2;
        KK = KK + J;
      }
    } else {
      JX = KX;
      JY = KY;
      for (J = 1; J <= N; J++) {
        TEMP1 = ALPHA * X[JX];
        TEMP2 = ZERO;
        IX = KX;
        IY = KY;
        for (K = KK; K <= KK + J - 2; K++) {
          Y[IY] = Y[IY] + TEMP1 * AP[K];
          TEMP2 = TEMP2 + AP[K] * X[IX];
          IX = IX + INCX;
          IY = IY + INCY;
        }
        Y[JY] = Y[JY] + TEMP1 * AP[KK + J - 1] + ALPHA * TEMP2;
        JX = JX + INCX;
        JY = JY + INCY;
        KK = KK + J;
      }
    }
  } else {
    // Form  y  when AP contains the lower triangle.

    if ((INCX == 1) && (INCY == 1)) {
      for (J = 1; J <= N; J++) {
        TEMP1 = ALPHA * X[J];
        TEMP2 = ZERO;
        Y[J] = Y[J] + TEMP1 * AP[KK];
        K = KK + 1;
        for (I = J + 1; I <= N; I++) {
          Y[I] = Y[I] + TEMP1 * AP[K];
          TEMP2 = TEMP2 + AP[K] * X[I];
          K = K + 1;
        }
        Y[J] = Y[J] + ALPHA * TEMP2;
        KK = KK + (N - J + 1);
      }
    } else {
      JX = KX;
      JY = KY;
      for (J = 1; J <= N; J++) {
        TEMP1 = ALPHA * X[JX];
        TEMP2 = ZERO;
        Y[JY] = Y[JY] + TEMP1 * AP[KK];
        IX = JX;
        IY = JY;
        for (K = KK + 1; K <= KK + N - J; K++) {
          IX = IX + INCX;
          IY = IY + INCY;
          Y[IY] = Y[IY] + TEMP1 * AP[K];
          TEMP2 = TEMP2 + AP[K] * X[IX];
        }
        Y[JY] = Y[JY] + ALPHA * TEMP2;
        JX = JX + INCX;
        JY = JY + INCY;
        KK = KK + (N - J + 1);
      }
    }
  }
}
