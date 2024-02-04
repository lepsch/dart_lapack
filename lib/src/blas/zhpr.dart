import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zhpr(
  final String UPLO,
  final int N,
  final double ALPHA,
  final Array<Complex> X,
  final int INCX,
  final Array<Complex> AP,
) {
// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  Complex TEMP;
  int I, INFO, IX, J, JX, K, KK = 0, KX = 0;

  // Test the input parameters.

  INFO = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO = 1;
  } else if (N < 0) {
    INFO = 2;
  } else if (INCX == 0) {
    INFO = 5;
  }
  if (INFO != 0) {
    xerbla('ZHPR  ', INFO);
    return;
  }

  // Quick return if possible.

  if ((N == 0) || (ALPHA == Complex.zero.toDouble())) return;

  // Set the start point in X if the increment is not unity.

  if (INCX <= 0) {
    KX = 1 - (N - 1) * INCX;
  } else if (INCX != 1) {
    KX = 1;
  }

  // Start the operations. In this version the elements of the array AP
  // are accessed sequentially with one pass through AP.

  KK = 1;
  if (lsame(UPLO, 'U')) {
    // Form  A  when upper triangle is stored in AP.

    if (INCX == 1) {
      for (J = 1; J <= N; J++) {
        if (X[J] != Complex.zero) {
          TEMP = ALPHA.toComplex() * X[J].conjugate();
          K = KK;
          for (I = 1; I <= J - 1; I++) {
            AP[K] = AP[K] + X[I] * TEMP;
            K = K + 1;
          }
          AP[KK + J - 1] =
              (AP[KK + J - 1].toDouble() + (X[J] * TEMP).toDouble())
                  .toComplex();
        } else {
          AP[KK + J - 1] = AP[KK + J - 1].real.toComplex();
        }
        KK = KK + J;
      }
    } else {
      JX = KX;
      for (J = 1; J <= N; J++) {
        if (X[JX] != Complex.zero) {
          TEMP = ALPHA.toComplex() * X[JX].conjugate();
          IX = KX;
          for (K = KK; K <= KK + J - 2; K++) {
            AP[K] = AP[K] + X[IX] * TEMP;
            IX = IX + INCX;
          }
          AP[KK + J - 1] =
              (AP[KK + J - 1].toDouble() + (X[JX] * TEMP).toDouble())
                  .toComplex();
        } else {
          AP[KK + J - 1] = (AP[KK + J - 1]).real.toComplex();
        }
        JX = JX + INCX;
        KK = KK + J;
      }
    }
  } else {
    // Form  A  when lower triangle is stored in AP.

    if (INCX == 1) {
      for (J = 1; J <= N; J++) {
        if (X[J] != Complex.zero) {
          TEMP = ALPHA.toComplex() * X[J].conjugate();
          AP[KK] = (AP[KK].toDouble() + (TEMP * X[J]).toDouble()).toComplex();
          K = KK + 1;
          for (I = J + 1; I <= N; I++) {
            AP[K] = AP[K] + X[I] * TEMP;
            K = K + 1;
          }
        } else {
          AP[KK] = AP[KK].real.toComplex();
        }
        KK = KK + N - J + 1;
      }
    } else {
      JX = KX;
      for (J = 1; J <= N; J++) {
        if (X[JX] != Complex.zero) {
          TEMP = ALPHA.toComplex() * X[JX].conjugate();
          AP[KK] = (AP[KK].toDouble() + (TEMP * X[JX]).toDouble()).toComplex();
          IX = JX;
          for (K = KK + 1; K <= KK + N - J; K++) {
            IX = IX + INCX;
            AP[K] = AP[K] + X[IX] * TEMP;
          }
        } else {
          AP[KK] = AP[KK].real.toComplex();
        }
        JX = JX + INCX;
        KK = KK + N - J + 1;
      }
    }
  }

  return;
}