// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';

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
  int INFO, IX, IY, JX, JY, KPLUS1, KX, KY, L;
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
    xerbla('ZSBMV', INFO);
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
        for (var I = 1; I <= N; I++) {
          Y[I] = Complex.zero;
        }
      } else {
        for (var I = 1; I <= N; I++) {
          Y[I] = BETA * Y[I];
        }
      }
    } else {
      IY = KY;
      if (BETA == Complex.zero) {
        for (var I = 1; I <= N; I++) {
          Y[IY] = Complex.zero;
          IY += INCY;
        }
      } else {
        for (var I = 1; I <= N; I++) {
          Y[IY] = BETA * Y[IY];
          IY += INCY;
        }
      }
    }
  }
  if (ALPHA == Complex.zero) return;
  if (lsame(UPLO, 'U')) {
    // Form  y  when upper triangle of A is stored.

    KPLUS1 = K + 1;
    if ((INCX == 1) && (INCY == 1)) {
      for (var J = 1; J <= N; J++) {
        TEMP1 = ALPHA * X[J];
        TEMP2 = Complex.zero;
        L = KPLUS1 - J;
        for (var I = max(1, J - K); I <= J - 1; I++) {
          Y[I] += TEMP1 * A[L + I][J];
          TEMP2 += A[L + I][J] * X[I];
        }
        Y[J] += TEMP1 * A[KPLUS1][J] + ALPHA * TEMP2;
      }
    } else {
      JX = KX;
      JY = KY;
      for (var J = 1; J <= N; J++) {
        TEMP1 = ALPHA * X[JX];
        TEMP2 = Complex.zero;
        IX = KX;
        IY = KY;
        L = KPLUS1 - J;
        for (var I = max(1, J - K); I <= J - 1; I++) {
          Y[IY] += TEMP1 * A[L + I][J];
          TEMP2 += A[L + I][J] * X[IX];
          IX += INCX;
          IY += INCY;
        }
        Y[JY] += TEMP1 * A[KPLUS1][J] + ALPHA * TEMP2;
        JX += INCX;
        JY += INCY;
        if (J > K) {
          KX += INCX;
          KY += INCY;
        }
      }
    }
  } else {
    // Form  y  when lower triangle of A is stored.

    if ((INCX == 1) && (INCY == 1)) {
      for (var J = 1; J <= N; J++) {
        TEMP1 = ALPHA * X[J];
        TEMP2 = Complex.zero;
        Y[J] += TEMP1 * A[1][J];
        L = 1 - J;
        for (var I = J + 1; I <= min(N, J + K); I++) {
          Y[I] += TEMP1 * A[L + I][J];
          TEMP2 += A[L + I][J] * X[I];
        }
        Y[J] += ALPHA * TEMP2;
      }
    } else {
      JX = KX;
      JY = KY;
      for (var J = 1; J <= N; J++) {
        TEMP1 = ALPHA * X[JX];
        TEMP2 = Complex.zero;
        Y[JY] += TEMP1 * A[1][J];
        L = 1 - J;
        IX = JX;
        IY = JY;
        for (var I = J + 1; I <= min(N, J + K); I++) {
          IX += INCX;
          IY += INCY;
          Y[IY] += TEMP1 * A[L + I][J];
          TEMP2 += A[L + I][J] * X[IX];
        }
        Y[JY] += ALPHA * TEMP2;
        JX += INCX;
        JY += INCY;
      }
    }
  }
}
