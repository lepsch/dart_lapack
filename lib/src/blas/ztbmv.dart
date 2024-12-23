// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/lsame.dart';
import 'package:dart_lapack/src/blas/xerbla.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

void ztbmv(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int N,
  final int K,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> X_,
  final int INCX,
) {
  final A = A_.having(ld: LDA);
  final X = X_.having();
  Complex TEMP;
  int I, INFO, IX, J, JX, KPLUS1, KX = 0, L;
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
  } else if (K < 0) {
    INFO = 5;
  } else if (LDA < (K + 1)) {
    INFO = 7;
  } else if (INCX == 0) {
    INFO = 9;
  }
  if (INFO != 0) {
    xerbla('ZTBMV', INFO);
    return;
  }

  // Quick return if possible.

  if (N == 0) return;

  NOCONJ = lsame(TRANS, 'T');
  NOUNIT = lsame(DIAG, 'N');

  // Set up the start point in X if the increment is not unity. This
  // will be  ( N - 1 )*INCX   too small for descending loops.

  if (INCX <= 0) {
    KX = 1 - (N - 1) * INCX;
  } else if (INCX != 1) {
    KX = 1;
  }

  // Start the operations. In this version the elements of A are
  // accessed sequentially with one pass through A.

  if (lsame(TRANS, 'N')) {
    // Form  x := A*x.

    if (lsame(UPLO, 'U')) {
      KPLUS1 = K + 1;
      if (INCX == 1) {
        for (J = 1; J <= N; J++) {
          if (X[J] != Complex.zero) {
            TEMP = X[J];
            L = KPLUS1 - J;
            for (I = max(1, J - K); I <= J - 1; I++) {
              X[I] += TEMP * A[L + I][J];
            }
            if (NOUNIT) X[J] *= A[KPLUS1][J];
          }
        }
      } else {
        JX = KX;
        for (J = 1; J <= N; J++) {
          if (X[JX] != Complex.zero) {
            TEMP = X[JX];
            IX = KX;
            L = KPLUS1 - J;
            for (I = max(1, J - K); I <= J - 1; I++) {
              X[IX] += TEMP * A[L + I][J];
              IX += INCX;
            }
            if (NOUNIT) X[JX] *= A[KPLUS1][J];
          }
          JX += INCX;
          if (J > K) KX += INCX;
        }
      }
    } else {
      if (INCX == 1) {
        for (J = N; J >= 1; J--) {
          if (X[J] != Complex.zero) {
            TEMP = X[J];
            L = 1 - J;
            for (I = min(N, J + K); I >= J + 1; I--) {
              X[I] += TEMP * A[L + I][J];
            }
            if (NOUNIT) X[J] *= A[1][J];
          }
        }
      } else {
        KX += (N - 1) * INCX;
        JX = KX;
        for (J = N; J >= 1; J--) {
          if (X[JX] != Complex.zero) {
            TEMP = X[JX];
            IX = KX;
            L = 1 - J;
            for (I = min(N, J + K); I >= J + 1; I--) {
              X[IX] += TEMP * A[L + I][J];
              IX -= INCX;
            }
            if (NOUNIT) X[JX] *= A[1][J];
          }
          JX -= INCX;
          if ((N - J) >= K) KX -= INCX;
        }
      }
    }
  } else {
    // Form  x := A**T*x  or  x := A**H*x.

    if (lsame(UPLO, 'U')) {
      KPLUS1 = K + 1;
      if (INCX == 1) {
        for (J = N; J >= 1; J--) {
          TEMP = X[J];
          L = KPLUS1 - J;
          if (NOCONJ) {
            if (NOUNIT) TEMP *= A[KPLUS1][J];
            for (I = J - 1; I >= max(1, J - K); I--) {
              TEMP += A[L + I][J] * X[I];
            }
          } else {
            if (NOUNIT) TEMP *= A[KPLUS1][J].conjugate();
            for (I = J - 1; I >= max(1, J - K); I--) {
              TEMP += A[L + I][J].conjugate() * X[I];
            }
          }
          X[J] = TEMP;
        }
      } else {
        KX += (N - 1) * INCX;
        JX = KX;
        for (J = N; J >= 1; J--) {
          TEMP = X[JX];
          KX -= INCX;
          IX = KX;
          L = KPLUS1 - J;
          if (NOCONJ) {
            if (NOUNIT) TEMP *= A[KPLUS1][J];
            for (I = J - 1; I >= max(1, J - K); I--) {
              TEMP += A[L + I][J] * X[IX];
              IX -= INCX;
            }
          } else {
            if (NOUNIT) TEMP *= A[KPLUS1][J].conjugate();
            for (I = J - 1; I >= max(1, J - K); I--) {
              TEMP += A[L + I][J].conjugate() * X[IX];
              IX -= INCX;
            }
          }
          X[JX] = TEMP;
          JX -= INCX;
        }
      }
    } else {
      if (INCX == 1) {
        for (J = 1; J <= N; J++) {
          TEMP = X[J];
          L = 1 - J;
          if (NOCONJ) {
            if (NOUNIT) TEMP *= A[1][J];
            for (I = J + 1; I <= min(N, J + K); I++) {
              TEMP += A[L + I][J] * X[I];
            }
          } else {
            if (NOUNIT) TEMP *= A[1][J].conjugate();
            for (I = J + 1; I <= min(N, J + K); I++) {
              TEMP += A[L + I][J].conjugate() * X[I];
            }
          }
          X[J] = TEMP;
        }
      } else {
        JX = KX;
        for (J = 1; J <= N; J++) {
          TEMP = X[JX];
          KX += INCX;
          IX = KX;
          L = 1 - J;
          if (NOCONJ) {
            if (NOUNIT) TEMP *= A[1][J];
            for (I = J + 1; I <= min(N, J + K); I++) {
              TEMP += A[L + I][J] * X[IX];
              IX += INCX;
            }
          } else {
            if (NOUNIT) TEMP *= A[1][J].conjugate();
            for (I = J + 1; I <= min(N, J + K); I++) {
              TEMP += A[L + I][J].conjugate() * X[IX];
              IX += INCX;
            }
          }
          X[JX] = TEMP;
          JX += INCX;
        }
      }
    }
  }
}
