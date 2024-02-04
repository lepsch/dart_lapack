import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zhpmv(
  final String UPLO,
  final int N,
  final Complex ALPHA,
  final Array<Complex> AP,
  final Array<Complex> X,
  final int INCX,
  final Complex BETA,
  final Array<Complex> Y,
  final int INCY,
) {
// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  Complex TEMP1, TEMP2;
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
    xerbla('ZHPMV ', INFO);
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
          Y[I] = Y[I] + TEMP1 * AP[K];
          TEMP2 = TEMP2 + AP[K].conjugate() * X[I];
          K = K + 1;
        }
        Y[J] = Y[J] + TEMP1 * AP[KK + J - 1].real.toComplex() + ALPHA * TEMP2;
        KK = KK + J;
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
          Y[IY] = Y[IY] + TEMP1 * AP[K];
          TEMP2 = TEMP2 + AP[K].conjugate() * X[IX];
          IX = IX + INCX;
          IY = IY + INCY;
        }
        Y[JY] =
            Y[JY] + TEMP1 * AP[KK + J - 1].real.toComplex() + ALPHA * TEMP2;
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
        TEMP2 = Complex.zero;
        Y[J] = Y[J] + TEMP1 * AP[KK].real.toComplex();
        K = KK + 1;
        for (I = J + 1; I <= N; I++) {
          Y[I] = Y[I] + TEMP1 * AP[K];
          TEMP2 = TEMP2 + AP[K].conjugate() * X[I];
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
        TEMP2 = Complex.zero;
        Y[JY] = Y[JY] + TEMP1 * AP[KK].real.toComplex();
        IX = JX;
        IY = JY;
        for (K = KK + 1; K <= KK + N - J; K++) {
          IX = IX + INCX;
          IY = IY + INCY;
          Y[IY] = Y[IY] + TEMP1 * AP[K];
          TEMP2 = TEMP2 + AP[K].conjugate() * X[IX];
        }
        Y[JY] = Y[JY] + ALPHA * TEMP2;
        JX = JX + INCX;
        JY = JY + INCY;
        KK = KK + (N - J + 1);
      }
    }
  }

  return;
}
