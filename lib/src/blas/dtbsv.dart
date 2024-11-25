// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/lsame.dart';
import 'package:dart_lapack/src/blas/xerbla.dart';
import 'package:dart_lapack/src/matrix.dart';

void dtbsv(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int N,
  final int K,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> X_,
  final int INCX,
) {
  final A = A_.having(ld: LDA);
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
  } else if (K < 0) {
    INFO = 5;
  } else if (LDA < (K + 1)) {
    INFO = 7;
  } else if (INCX == 0) {
    INFO = 9;
  }
  if (INFO != 0) {
    xerbla('DTBSV', INFO);
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

  // Start the operations. In this version the elements of A are
  // accessed by sequentially with one pass through A.

  if (lsame(TRANS, 'N')) {
    // Form  x := inv( A )*x.

    if (lsame(UPLO, 'U')) {
      final KPLUS1 = K + 1;
      if (INCX == 1) {
        for (var J = N; J >= 1; J--) {
          if (X[J] != ZERO) {
            final L = KPLUS1 - J;
            if (NOUNIT) X[J] /= A[KPLUS1][J];
            final TEMP = X[J];
            for (var I = J - 1; I >= max(1, J - K); I--) {
              X[I] -= TEMP * A[L + I][J];
            }
          }
        }
      } else {
        KX += (N - 1) * INCX;
        var JX = KX;
        for (var J = N; J >= 1; J--) {
          KX -= INCX;
          if (X[JX] != ZERO) {
            var IX = KX;
            final L = KPLUS1 - J;
            if (NOUNIT) X[JX] /= A[KPLUS1][J];
            final TEMP = X[JX];
            for (var I = J - 1; I >= max(1, J - K); I--) {
              X[IX] -= TEMP * A[L + I][J];
              IX -= INCX;
            }
          }
          JX -= INCX;
        }
      }
    } else {
      if (INCX == 1) {
        for (var J = 1; J <= N; J++) {
          if (X[J] != ZERO) {
            final L = 1 - J;
            if (NOUNIT) X[J] /= A[1][J];
            final TEMP = X[J];
            for (var I = J + 1; I <= min(N, J + K); I++) {
              X[I] -= TEMP * A[L + I][J];
            }
          }
        }
      } else {
        var JX = KX;
        for (var J = 1; J <= N; J++) {
          KX += INCX;
          if (X[JX] != ZERO) {
            var IX = KX;
            final L = 1 - J;
            if (NOUNIT) X[JX] /= A[1][J];
            final TEMP = X[JX];
            for (var I = J + 1; I <= min(N, J + K); I++) {
              X[IX] -= TEMP * A[L + I][J];
              IX += INCX;
            }
          }
          JX += INCX;
        }
      }
    }
  } else {
    // Form  x := inv( A**T)*x.

    if (lsame(UPLO, 'U')) {
      final KPLUS1 = K + 1;
      if (INCX == 1) {
        for (var J = 1; J <= N; J++) {
          var TEMP = X[J];
          final L = KPLUS1 - J;
          for (var I = max(1, J - K); I <= J - 1; I++) {
            TEMP -= A[L + I][J] * X[I];
          }
          if (NOUNIT) TEMP /= A[KPLUS1][J];
          X[J] = TEMP;
        }
      } else {
        var JX = KX;
        for (var J = 1; J <= N; J++) {
          var TEMP = X[JX];
          var IX = KX;
          final L = KPLUS1 - J;
          for (var I = max(1, J - K); I <= J - 1; I++) {
            TEMP -= A[L + I][J] * X[IX];
            IX += INCX;
          }
          if (NOUNIT) TEMP /= A[KPLUS1][J];
          X[JX] = TEMP;
          JX += INCX;
          if (J > K) KX += INCX;
        }
      }
    } else {
      if (INCX == 1) {
        for (var J = N; J >= 1; J--) {
          var TEMP = X[J];
          final L = 1 - J;
          for (var I = min(N, J + K); I >= J + 1; I--) {
            TEMP -= A[L + I][J] * X[I];
          }
          if (NOUNIT) TEMP /= A[1][J];
          X[J] = TEMP;
        }
      } else {
        KX += (N - 1) * INCX;
        var JX = KX;
        for (var J = N; J >= 1; J--) {
          var TEMP = X[JX];
          var IX = KX;
          final L = 1 - J;
          for (var I = min(N, J + K); I >= J + 1; I--) {
            TEMP -= A[L + I][J] * X[IX];
            IX -= INCX;
          }
          if (NOUNIT) TEMP /= A[1][J];
          X[JX] = TEMP;
          JX -= INCX;
          if ((N - J) >= K) KX -= INCX;
        }
      }
    }
  }
}
