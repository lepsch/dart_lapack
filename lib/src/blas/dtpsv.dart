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

  // Test the input parameters.

  var INFO = 0;
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

  final NOUNIT = lsame(DIAG, 'N');

  // Set up the start point in X if the increment is not unity. This
  // will be  ( N - 1 )*INCX  too small for descending loops.

  var KX = switch (INCX) {
    <= 0 => 1 - (N - 1) * INCX,
    1 => 0,
    _ => 1,
  };

  // Start the operations. In this version the elements of AP are
  // accessed sequentially with one pass through AP.

  if (lsame(TRANS, 'N')) {
    // Form  x := inv( A )*x.

    if (lsame(UPLO, 'U')) {
      var KK = (N * (N + 1)) ~/ 2;
      if (INCX == 1) {
        for (var J = N; J >= 1; J--) {
          if (X[J] != ZERO) {
            if (NOUNIT) X[J] /= AP[KK];
            final TEMP = X[J];
            for (var I = J - 1, K = KK - 1; I >= 1; I--, K--) {
              X[I] -= TEMP * AP[K];
            }
          }
          KK -= J;
        }
      } else {
        var JX = KX + (N - 1) * INCX;
        for (var J = N; J >= 1; J--) {
          if (X[JX] != ZERO) {
            if (NOUNIT) X[JX] /= AP[KK];
            final TEMP = X[JX];
            var IX = JX;
            for (var K = KK - 1; K >= KK - J + 1; K--) {
              IX -= INCX;
              X[IX] -= TEMP * AP[K];
            }
          }
          JX -= INCX;
          KK -= J;
        }
      }
    } else {
      var KK = 1;
      if (INCX == 1) {
        for (var J = 1; J <= N; J++) {
          if (X[J] != ZERO) {
            if (NOUNIT) X[J] /= AP[KK];
            final TEMP = X[J];
            for (var I = J + 1, K = KK + 1; I <= N; I++, K++) {
              X[I] -= TEMP * AP[K];
            }
          }
          KK += N - J + 1;
        }
      } else {
        var JX = KX;
        for (var J = 1; J <= N; J++) {
          if (X[JX] != ZERO) {
            if (NOUNIT) X[JX] /= AP[KK];
            final TEMP = X[JX];
            var IX = JX;
            for (var K = KK + 1; K <= KK + N - J; K++) {
              IX += INCX;
              X[IX] -= TEMP * AP[K];
            }
          }
          JX += INCX;
          KK += N - J + 1;
        }
      }
    }
  } else {
    // Form  x := inv( A**T )*x.

    if (lsame(UPLO, 'U')) {
      var KK = 1;
      if (INCX == 1) {
        for (var J = 1; J <= N; J++) {
          var TEMP = X[J];
          for (var I = 1, K = KK; I <= J - 1; I++, K++) {
            TEMP -= AP[K] * X[I];
          }
          if (NOUNIT) TEMP /= AP[KK + J - 1];
          X[J] = TEMP;
          KK += J;
        }
      } else {
        var JX = KX;
        for (var J = 1; J <= N; J++) {
          var TEMP = X[JX];
          var IX = KX;
          for (var K = KK; K <= KK + J - 2; K++) {
            TEMP -= AP[K] * X[IX];
            IX += INCX;
          }
          if (NOUNIT) TEMP /= AP[KK + J - 1];
          X[JX] = TEMP;
          JX += INCX;
          KK += J;
        }
      }
    } else {
      var KK = (N * (N + 1)) ~/ 2;
      if (INCX == 1) {
        for (var J = N; J >= 1; J--) {
          var TEMP = X[J];
          for (var I = N, K = KK; I >= J + 1; I--, K--) {
            TEMP -= AP[K] * X[I];
          }
          if (NOUNIT) TEMP /= AP[KK - N + J];
          X[J] = TEMP;
          KK -= N - J + 1;
        }
      } else {
        KX += (N - 1) * INCX;
        var JX = KX;
        for (var J = N; J >= 1; J--) {
          var TEMP = X[JX];
          var IX = KX;
          for (var K = KK; K >= KK - (N - (J + 1)); K--) {
            TEMP -= AP[K] * X[IX];
            IX -= INCX;
          }
          if (NOUNIT) TEMP /= AP[KK - N + J];
          X[JX] = TEMP;
          JX -= INCX;
          KK -= N - J + 1;
        }
      }
    }
  }
}
