import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zspmv(
  final String UPLO,
  final int N,
  final Complex ALPHA,
  final Array<Complex> AP_,
  final Array<Complex> X_,
  final int INCX,
  final Complex BETA,
  final Array<Complex> Y_,
  final int INCY,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.having();
  final X = X_.having();
  final Y = Y_.having();
  int I, INFO, IX, IY, J, JX, JY, K, KK, KX, KY;
  Complex TEMP1, TEMP2;

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
    xerbla('ZSPMV', INFO);
    return;
  }

  // Quick return if possible.

  if ((N == 0) || ((ALPHA == Complex.zero) && (BETA == Complex.one))) return;

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

  if (BETA != Complex.one) {
    if (INCY == 1) {
      if (BETA == Complex.zero) {
        for (I = 1; I <= N; I++) {
          Y[I] = Complex.zero;
        }
      } else {
        for (I = 1; I <= N; I++) {
          Y[I] = BETA * Y[I];
        }
      }
    } else {
      IY = KY;
      if (BETA == Complex.zero) {
        for (I = 1; I <= N; I++) {
          Y[IY] = Complex.zero;
          IY += INCY;
        }
      } else {
        for (I = 1; I <= N; I++) {
          Y[IY] = BETA * Y[IY];
          IY += INCY;
        }
      }
    }
  }
  if (ALPHA == Complex.zero) return;
  KK = 1;
  if (lsame(UPLO, 'U')) {
    // Form  y  when AP contains the upper triangle.

    if ((INCX == 1) && (INCY == 1)) {
      for (J = 1; J <= N; J++) {
        TEMP1 = ALPHA * X[J];
        TEMP2 = Complex.zero;
        K = KK;
        for (I = 1; I <= J - 1; I++) {
          Y[I] += TEMP1 * AP[K];
          TEMP2 += AP[K] * X[I];
          K++;
        }
        Y[J] += TEMP1 * AP[KK + J - 1] + ALPHA * TEMP2;
        KK += J;
      }
    } else {
      JX = KX;
      JY = KY;
      for (J = 1; J <= N; J++) {
        TEMP1 = ALPHA * X[JX];
        TEMP2 = Complex.zero;
        IX = KX;
        IY = KY;
        for (K = KK; K <= KK + J - 2; K++) {
          Y[IY] += TEMP1 * AP[K];
          TEMP2 += AP[K] * X[IX];
          IX += INCX;
          IY += INCY;
        }
        Y[JY] += TEMP1 * AP[KK + J - 1] + ALPHA * TEMP2;
        JX += INCX;
        JY += INCY;
        KK += J;
      }
    }
  } else {
    // Form  y  when AP contains the lower triangle.

    if ((INCX == 1) && (INCY == 1)) {
      for (J = 1; J <= N; J++) {
        TEMP1 = ALPHA * X[J];
        TEMP2 = Complex.zero;
        Y[J] += TEMP1 * AP[KK];
        K = KK + 1;
        for (I = J + 1; I <= N; I++) {
          Y[I] += TEMP1 * AP[K];
          TEMP2 += AP[K] * X[I];
          K++;
        }
        Y[J] += ALPHA * TEMP2;
        KK += (N - J + 1);
      }
    } else {
      JX = KX;
      JY = KY;
      for (J = 1; J <= N; J++) {
        TEMP1 = ALPHA * X[JX];
        TEMP2 = Complex.zero;
        Y[JY] += TEMP1 * AP[KK];
        IX = JX;
        IY = JY;
        for (K = KK + 1; K <= KK + N - J; K++) {
          IX += INCX;
          IY += INCY;
          Y[IY] += TEMP1 * AP[K];
          TEMP2 += AP[K] * X[IX];
        }
        Y[JY] += ALPHA * TEMP2;
        JX += INCX;
        JY += INCY;
        KK += (N - J + 1);
      }
    }
  }
}
