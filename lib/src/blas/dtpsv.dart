import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/matrix.dart';

void dtpsv(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int N,
  final Array<double> AP_,
  final Array<double> X_,
  final int INCX,
) {
// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.having();
  final X = X_.having();
  const ZERO = 0.0;
  double TEMP;
  int I, INFO, IX, J, JX, K, KK, KX = 0;
  bool NOUNIT;

  // Test the input parameters.

  INFO = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO = 1;
  } else if (!lsame(TRANS, 'N') && !lsame(TRANS, 'T') && !lsame(TRANS, 'C')) {
    INFO = 2;
  } else if (!lsame(DIAG, 'U') && !lsame(DIAG, 'N')) {
    INFO = 3;
  } else if (N < 0) {
    INFO = 4;
  } else if (INCX == 0) {
    INFO = 7;
  }
  if (INFO != 0) {
    xerbla('DTPSV ', INFO);
    return;
  }

  // Quick return if possible.

  if (N == 0) return;

  NOUNIT = lsame(DIAG, 'N');

  // Set up the start point in X if the increment is not unity. This
  // will be  ( N - 1 )*INCX  too small for descending loops.

  if (INCX <= 0) {
    KX = 1 - (N - 1) * INCX;
  } else if (INCX != 1) {
    KX = 1;
  }

  // Start the operations. In this version the elements of AP are
  // accessed sequentially with one pass through AP.

  if (lsame(TRANS, 'N')) {
    // Form  x := inv( A )*x.

    if (lsame(UPLO, 'U')) {
      KK = (N * (N + 1)) ~/ 2;
      if (INCX == 1) {
        for (J = N; J >= 1; J--) {
          if (X[J] != ZERO) {
            if (NOUNIT) X[J] = X[J] / AP[KK];
            TEMP = X[J];
            K = KK - 1;
            for (I = J - 1; I >= 1; I--) {
              X[I] = X[I] - TEMP * AP[K];
              K = K - 1;
            }
          }
          KK = KK - J;
        }
      } else {
        JX = KX + (N - 1) * INCX;
        for (J = N; J >= 1; J--) {
          if (X[JX] != ZERO) {
            if (NOUNIT) X[JX] = X[JX] / AP[KK];
            TEMP = X[JX];
            IX = JX;
            for (K = KK - 1; K >= KK - J + 1; K--) {
              IX = IX - INCX;
              X[IX] = X[IX] - TEMP * AP[K];
            }
          }
          JX = JX - INCX;
          KK = KK - J;
        }
      }
    } else {
      KK = 1;
      if (INCX == 1) {
        for (J = 1; J <= N; J++) {
          if (X[J] != ZERO) {
            if (NOUNIT) X[J] = X[J] / AP[KK];
            TEMP = X[J];
            K = KK + 1;
            for (I = J + 1; I <= N; I++) {
              X[I] = X[I] - TEMP * AP[K];
              K = K + 1;
            }
          }
          KK = KK + (N - J + 1);
        }
      } else {
        JX = KX;
        for (J = 1; J <= N; J++) {
          if (X[JX] != ZERO) {
            if (NOUNIT) X[JX] = X[JX] / AP[KK];
            TEMP = X[JX];
            IX = JX;
            for (K = KK + 1; K <= KK + N - J; K++) {
              IX = IX + INCX;
              X[IX] = X[IX] - TEMP * AP[K];
            }
          }
          JX = JX + INCX;
          KK = KK + (N - J + 1);
        }
      }
    }
  } else {
    // Form  x := inv( A**T )*x.

    if (lsame(UPLO, 'U')) {
      KK = 1;
      if (INCX == 1) {
        for (J = 1; J <= N; J++) {
          TEMP = X[J];
          K = KK;
          for (I = 1; I <= J - 1; I++) {
            TEMP = TEMP - AP[K] * X[I];
            K = K + 1;
          }
          if (NOUNIT) TEMP = TEMP / AP[KK + J - 1];
          X[J] = TEMP;
          KK = KK + J;
        }
      } else {
        JX = KX;
        for (J = 1; J <= N; J++) {
          TEMP = X[JX];
          IX = KX;
          for (K = KK; K <= KK + J - 2; K++) {
            TEMP = TEMP - AP[K] * X[IX];
            IX = IX + INCX;
          }
          if (NOUNIT) TEMP = TEMP / AP[KK + J - 1];
          X[JX] = TEMP;
          JX = JX + INCX;
          KK = KK + J;
        }
      }
    } else {
      KK = (N * (N + 1)) ~/ 2;
      if (INCX == 1) {
        for (J = N; J >= 1; J--) {
          TEMP = X[J];
          K = KK;
          for (I = N; I >= J + 1; I--) {
            TEMP = TEMP - AP[K] * X[I];
            K = K - 1;
          }
          if (NOUNIT) TEMP = TEMP / AP[KK - N + J];
          X[J] = TEMP;
          KK = KK - (N - J + 1);
        }
      } else {
        KX = KX + (N - 1) * INCX;
        JX = KX;
        for (J = N; J >= 1; J--) {
          TEMP = X[JX];
          IX = KX;
          for (K = KK; K >= KK - (N - (J + 1)); K--) {
            TEMP = TEMP - AP[K] * X[IX];
            IX = IX - INCX;
          }
          if (NOUNIT) TEMP = TEMP / AP[KK - N + J];
          X[JX] = TEMP;
          JX = JX - INCX;
          KK = KK - (N - J + 1);
        }
      }
    }
  }
}
