import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zhpr2(
  final String UPLO,
  final int N,
  final Complex ALPHA,
  final Array<Complex> X,
  final int INCX,
  final Array<Complex> Y,
  final int INCY,
  final Array<Complex> AP,
) {
// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  Complex TEMP1, TEMP2;
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
    xerbla('ZHPR2 ', INFO);
    return;
  }

  // Quick return if possible.

  if ((N == 0) || (ALPHA == Complex.zero)) return;

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
        if ((X[J] != Complex.zero) || (Y[J] != Complex.zero)) {
          TEMP1 = ALPHA * Y[J].conjugate();
          TEMP2 = ALPHA * X[J].conjugate();
          K = KK;
          for (I = 1; I <= J - 1; I++) {
            AP[K] = AP[K] + X[I] * TEMP1 + Y[I] * TEMP2;
            K = K + 1;
          }
          AP[KK + J - 1] = (AP[KK + J - 1].toDouble() +
                  (X[J] * TEMP1 + Y[J] * TEMP2).toDouble())
              .toComplex();
        } else {
          AP[KK + J - 1] = AP[KK + J - 1].real.toComplex();
        }
        KK = KK + J;
      }
    } else {
      for (J = 1; J <= N; J++) {
        if ((X[JX] != Complex.zero) || (Y[JY] != Complex.zero)) {
          TEMP1 = ALPHA * Y[JY].conjugate();
          TEMP2 = (ALPHA * X[JX]).conjugate();
          IX = KX;
          IY = KY;
          for (K = KK; K <= KK + J - 2; K++) {
            AP[K] = AP[K] + X[IX] * TEMP1 + Y[IY] * TEMP2;
            IX = IX + INCX;
            IY = IY + INCY;
          }
          AP[KK + J - 1] = (AP[KK + J - 1].toDouble() +
                  (X[JX] * TEMP1 + Y[JY] * TEMP2).toDouble())
              .toComplex();
        } else {
          AP[KK + J - 1] = AP[KK + J - 1].real.toComplex();
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
        if ((X[J] != Complex.zero) || (Y[J] != Complex.zero)) {
          TEMP1 = ALPHA * Y[J].conjugate();
          TEMP2 = ALPHA * X[J].conjugate();
          AP[KK] =
              (AP[KK].toDouble() + (X[J] * TEMP1 + Y[J] * TEMP2).toDouble())
                  .toComplex();
          K = KK + 1;
          for (I = J + 1; I <= N; I++) {
            AP[K] = AP[K] + X[I] * TEMP1 + Y[I] * TEMP2;
            K = K + 1;
          }
        } else {
          AP[KK] = AP[KK].real.toComplex();
        }
        KK = KK + N - J + 1;
      }
    } else {
      for (J = 1; J <= N; J++) {
        if ((X[JX] != Complex.zero) || (Y[JY] != Complex.zero)) {
          TEMP1 = ALPHA * Y[JY].conjugate();
          TEMP2 = (ALPHA * X[JX]).conjugate();
          AP[KK] =
              (AP[KK].toDouble() + (X[JX] * TEMP1 + Y[JY] * TEMP2).toDouble())
                  .toComplex();
          IX = JX;
          IY = JY;
          for (K = KK + 1; K <= KK + N - J; K++) {
            IX = IX + INCX;
            IY = IY + INCY;
            AP[K] = AP[K] + X[IX] * TEMP1 + Y[IY] * TEMP2;
          }
        } else {
          AP[KK] = AP[KK].real.toComplex();
        }
        JX = JX + INCX;
        JY = JY + INCY;
        KK = KK + N - J + 1;
      }
    }
  }

  return;
}
