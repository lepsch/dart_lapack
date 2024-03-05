import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/matrix.dart';

void dspr(
  final String UPLO,
  final int N,
  final double ALPHA,
  final Array<double> X_,
  final int INCX,
  final Array<double> AP,
) {
// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having();
  const ZERO = 0.0;
  double TEMP;
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
    xerbla('DSPR  ', INFO);
    return;
  }

  // Quick return if possible.

  if ((N == 0) || (ALPHA == ZERO)) return;

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
        if (X[J] != ZERO) {
          TEMP = ALPHA * X[J];
          K = KK;
          for (I = 1; I <= J; I++) {
            AP[K] = AP[K] + X[I] * TEMP;
            K = K + 1;
          }
        }
        KK = KK + J;
      }
    } else {
      JX = KX;
      for (J = 1; J <= N; J++) {
        if (X[JX] != ZERO) {
          TEMP = ALPHA * X[JX];
          IX = KX;
          for (K = KK; K <= KK + J - 1; K++) {
            AP[K] = AP[K] + X[IX] * TEMP;
            IX = IX + INCX;
          }
        }
        JX = JX + INCX;
        KK = KK + J;
      }
    }
  } else {
    // Form  A  when lower triangle is stored in AP.

    if (INCX == 1) {
      for (J = 1; J <= N; J++) {
        if (X[J] != ZERO) {
          TEMP = ALPHA * X[J];
          K = KK;
          for (I = J; I <= N; I++) {
            AP[K] = AP[K] + X[I] * TEMP;
            K = K + 1;
          }
        }
        KK = KK + N - J + 1;
      }
    } else {
      JX = KX;
      for (J = 1; J <= N; J++) {
        if (X[JX] != ZERO) {
          TEMP = ALPHA * X[JX];
          IX = JX;
          for (K = KK; K <= KK + N - J; K++) {
            AP[K] = AP[K] + X[IX] * TEMP;
            IX = IX + INCX;
          }
        }
        JX = JX + INCX;
        KK = KK + N - J + 1;
      }
    }
  }
}
