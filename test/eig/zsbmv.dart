import 'dart:math';

import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zsbmv(
  final String UPLO,
  final int N,
  final int K,
  final Complex ALPHA,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> X_,
  final int INCX,
  final Complex BETA,
  final Array<Complex> Y_,
  final int INCY,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final X = X_.having();
  final Y = Y_.having();
  int I, INFO, IX, IY, J, JX, JY, KPLUS1, KX, KY, L;
  Complex TEMP1, TEMP2;

  // Test the input parameters.

  INFO = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO = 1;
  } else if (N < 0) {
    INFO = 2;
  } else if (K < 0) {
    INFO = 3;
  } else if (LDA < (K + 1)) {
    INFO = 6;
  } else if (INCX == 0) {
    INFO = 8;
  } else if (INCY == 0) {
    INFO = 11;
  }
  if (INFO != 0) {
    xerbla('ZSBMV ', INFO);
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

  // Start the operations. In this version the elements of the array A
  // are accessed sequentially with one pass through A.

  // First form  y := beta*y.

  if (BETA != Complex.one) {
    if (INCY == 1) {
      if (BETA == Complex.zero) {
        for (I = 1; I <= N; I++) {
          // 10
          Y[I] = Complex.zero;
        } // 10
      } else {
        for (I = 1; I <= N; I++) {
          // 20
          Y[I] = BETA * Y[I];
        } // 20
      }
    } else {
      IY = KY;
      if (BETA == Complex.zero) {
        for (I = 1; I <= N; I++) {
          // 30
          Y[IY] = Complex.zero;
          IY += INCY;
        } // 30
      } else {
        for (I = 1; I <= N; I++) {
          // 40
          Y[IY] = BETA * Y[IY];
          IY += INCY;
        } // 40
      }
    }
  }
  if (ALPHA == Complex.zero) return;
  if (lsame(UPLO, 'U')) {
    // Form  y  when upper triangle of A is stored.

    KPLUS1 = K + 1;
    if ((INCX == 1) && (INCY == 1)) {
      for (J = 1; J <= N; J++) {
        // 60
        TEMP1 = ALPHA * X[J];
        TEMP2 = Complex.zero;
        L = KPLUS1 - J;
        for (I = max(1, J - K); I <= J - 1; I++) {
          // 50
          Y[I] = Y[I] + TEMP1 * A[L + I][J];
          TEMP2 += A[L + I][J] * X[I];
        } // 50
        Y[J] = Y[J] + TEMP1 * A[KPLUS1][J] + ALPHA * TEMP2;
      } // 60
    } else {
      JX = KX;
      JY = KY;
      for (J = 1; J <= N; J++) {
        // 80
        TEMP1 = ALPHA * X[JX];
        TEMP2 = Complex.zero;
        IX = KX;
        IY = KY;
        L = KPLUS1 - J;
        for (I = max(1, J - K); I <= J - 1; I++) {
          // 70
          Y[IY] = Y[IY] + TEMP1 * A[L + I][J];
          TEMP2 += A[L + I][J] * X[IX];
          IX += INCX;
          IY += INCY;
        } // 70
        Y[JY] = Y[JY] + TEMP1 * A[KPLUS1][J] + ALPHA * TEMP2;
        JX += INCX;
        JY += INCY;
        if (J > K) {
          KX += INCX;
          KY += INCY;
        }
      } // 80
    }
  } else {
    // Form  y  when lower triangle of A is stored.

    if ((INCX == 1) && (INCY == 1)) {
      for (J = 1; J <= N; J++) {
        // 100
        TEMP1 = ALPHA * X[J];
        TEMP2 = Complex.zero;
        Y[J] = Y[J] + TEMP1 * A[1][J];
        L = 1 - J;
        for (I = J + 1; I <= min(N, J + K); I++) {
          // 90
          Y[I] = Y[I] + TEMP1 * A[L + I][J];
          TEMP2 += A[L + I][J] * X[I];
        } // 90
        Y[J] = Y[J] + ALPHA * TEMP2;
      } // 100
    } else {
      JX = KX;
      JY = KY;
      for (J = 1; J <= N; J++) {
        // 120
        TEMP1 = ALPHA * X[JX];
        TEMP2 = Complex.zero;
        Y[JY] = Y[JY] + TEMP1 * A[1][J];
        L = 1 - J;
        IX = JX;
        IY = JY;
        for (I = J + 1; I <= min(N, J + K); I++) {
          // 110
          IX += INCX;
          IY += INCY;
          Y[IY] = Y[IY] + TEMP1 * A[L + I][J];
          TEMP2 += A[L + I][J] * X[IX];
        } // 110
        Y[JY] = Y[JY] + ALPHA * TEMP2;
        JX += INCX;
        JY += INCY;
      } // 120
    }
  }
}
