import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void ztpmv(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int N,
  final Array<Complex> AP_,
  final Array<Complex> X_,
  final int INCX,
) {
// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.having();
  final X = X_.having();
  Complex TEMP;
  int I, INFO, IX, J, JX, K, KK, KX = 0;
  bool NOCONJ, NOUNIT;

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
    xerbla('ZTPMV ', INFO);
    return;
  }

  // Quick return if possible.

  if (N == 0) return;

  NOCONJ = lsame(TRANS, 'T');
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
    // Form  x:= A*x.

    if (lsame(UPLO, 'U')) {
      KK = 1;
      if (INCX == 1) {
        for (J = 1; J <= N; J++) {
          if (X[J] != Complex.zero) {
            TEMP = X[J];
            K = KK;
            for (I = 1; I <= J - 1; I++) {
              X[I] += TEMP * AP[K];
              K++;
            }
            if (NOUNIT) X[J] *= AP[KK + J - 1];
          }
          KK += J;
        }
      } else {
        JX = KX;
        for (J = 1; J <= N; J++) {
          if (X[JX] != Complex.zero) {
            TEMP = X[JX];
            IX = KX;
            for (K = KK; K <= KK + J - 2; K++) {
              X[IX] += TEMP * AP[K];
              IX += INCX;
            }
            if (NOUNIT) X[JX] *= AP[KK + J - 1];
          }
          JX += INCX;
          KK += J;
        }
      }
    } else {
      KK = (N * (N + 1)) ~/ 2;
      if (INCX == 1) {
        for (J = N; J >= 1; J--) {
          if (X[J] != Complex.zero) {
            TEMP = X[J];
            K = KK;
            for (I = N; I >= J + 1; I--) {
              X[I] += TEMP * AP[K];
              K--;
            }
            if (NOUNIT) X[J] *= AP[KK - N + J];
          }
          KK -= (N - J + 1);
        }
      } else {
        KX += (N - 1) * INCX;
        JX = KX;
        for (J = N; J >= 1; J--) {
          if (X[JX] != Complex.zero) {
            TEMP = X[JX];
            IX = KX;
            for (K = KK; K >= KK - (N - (J + 1)); K--) {
              X[IX] += TEMP * AP[K];
              IX -= INCX;
            }
            if (NOUNIT) X[JX] *= AP[KK - N + J];
          }
          JX -= INCX;
          KK -= (N - J + 1);
        }
      }
    }
  } else {
    // Form  x := A**T*x  or  x := A**H*x.

    if (lsame(UPLO, 'U')) {
      KK = (N * (N + 1)) ~/ 2;
      if (INCX == 1) {
        for (J = N; J >= 1; J--) {
          TEMP = X[J];
          K = KK - 1;
          if (NOCONJ) {
            if (NOUNIT) TEMP = TEMP * AP[KK];
            for (I = J - 1; I >= 1; I--) {
              TEMP += AP[K] * X[I];
              K--;
            }
          } else {
            if (NOUNIT) TEMP = TEMP * AP[KK].conjugate();
            for (I = J - 1; I >= 1; I--) {
              TEMP += AP[K].conjugate() * X[I];
              K--;
            }
          }
          X[J] = TEMP;
          KK -= J;
        }
      } else {
        JX = KX + (N - 1) * INCX;
        for (J = N; J >= 1; J--) {
          TEMP = X[JX];
          IX = JX;
          if (NOCONJ) {
            if (NOUNIT) TEMP = TEMP * AP[KK];
            for (K = KK - 1; K >= KK - J + 1; K--) {
              IX -= INCX;
              TEMP += AP[K] * X[IX];
            }
          } else {
            if (NOUNIT) TEMP = TEMP * AP[KK].conjugate();
            for (K = KK - 1; K >= KK - J + 1; K--) {
              IX -= INCX;
              TEMP += AP[K].conjugate() * X[IX];
            }
          }
          X[JX] = TEMP;
          JX -= INCX;
          KK -= J;
        }
      }
    } else {
      KK = 1;
      if (INCX == 1) {
        for (J = 1; J <= N; J++) {
          TEMP = X[J];
          K = KK + 1;
          if (NOCONJ) {
            if (NOUNIT) TEMP = TEMP * AP[KK];
            for (I = J + 1; I <= N; I++) {
              TEMP += AP[K] * X[I];
              K++;
            }
          } else {
            if (NOUNIT) TEMP = TEMP * AP[KK].conjugate();
            for (I = J + 1; I <= N; I++) {
              TEMP += AP[K].conjugate() * X[I];
              K++;
            }
          }
          X[J] = TEMP;
          KK += (N - J + 1);
        }
      } else {
        JX = KX;
        for (J = 1; J <= N; J++) {
          TEMP = X[JX];
          IX = JX;
          if (NOCONJ) {
            if (NOUNIT) TEMP = TEMP * AP[KK];
            for (K = KK + 1; K <= KK + N - J; K++) {
              IX += INCX;
              TEMP += AP[K] * X[IX];
            }
          } else {
            if (NOUNIT) TEMP = TEMP * AP[KK].conjugate();
            for (K = KK + 1; K <= KK + N - J; K++) {
              IX += INCX;
              TEMP += AP[K].conjugate() * X[IX];
            }
          }
          X[JX] = TEMP;
          JX += INCX;
          KK += (N - J + 1);
        }
      }
    }
  }
}
