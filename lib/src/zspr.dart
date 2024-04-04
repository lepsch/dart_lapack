import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zspr(
  final String UPLO,
  final int N,
  final Complex ALPHA,
  final Array<Complex> X_,
  final int INCX,
  final Array<Complex> AP_,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having();
  final AP = AP_.having();
  int I, INFO, IX, J, JX, K, KK, KX = 0;
  Complex TEMP;

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
    xerbla('ZSPR', INFO);
    return;
  }

  // Quick return if possible.

  if ((N == 0) || (ALPHA == Complex.zero)) return;

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
          TEMP = ALPHA * X[J];
          K = KK;
          for (I = 1; I <= J - 1; I++) {
            AP[K] += X[I] * TEMP;
            K++;
          }
          AP[KK + J - 1] += X[J] * TEMP;
        } else {
          AP[KK + J - 1] = AP[KK + J - 1];
        }
        KK += J;
      }
    } else {
      JX = KX;
      for (J = 1; J <= N; J++) {
        if (X[JX] != Complex.zero) {
          TEMP = ALPHA * X[JX];
          IX = KX;
          for (K = KK; K <= KK + J - 2; K++) {
            AP[K] += X[IX] * TEMP;
            IX += INCX;
          }
          AP[KK + J - 1] += X[JX] * TEMP;
        } else {
          AP[KK + J - 1] = AP[KK + J - 1];
        }
        JX += INCX;
        KK += J;
      }
    }
  } else {
    // Form  A  when lower triangle is stored in AP.

    if (INCX == 1) {
      for (J = 1; J <= N; J++) {
        if (X[J] != Complex.zero) {
          TEMP = ALPHA * X[J];
          AP[KK] += TEMP * X[J];
          K = KK + 1;
          for (I = J + 1; I <= N; I++) {
            AP[K] += X[I] * TEMP;
            K++;
          }
        } else {
          AP[KK] = AP[KK];
        }
        KK += N - J + 1;
      }
    } else {
      JX = KX;
      for (J = 1; J <= N; J++) {
        if (X[JX] != Complex.zero) {
          TEMP = ALPHA * X[JX];
          AP[KK] += TEMP * X[JX];
          IX = JX;
          for (K = KK + 1; K <= KK + N - J; K++) {
            IX += INCX;
            AP[K] += X[IX] * TEMP;
          }
        } else {
          AP[KK] = AP[KK];
        }
        JX += INCX;
        KK += N - J + 1;
      }
    }
  }
}
