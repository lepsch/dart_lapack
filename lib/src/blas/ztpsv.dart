// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void ztpsv(
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
    xerbla('ZTPSV', INFO);
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
    // Form  x := inv( A )*x.

    if (lsame(UPLO, 'U')) {
      KK = (N * (N + 1)) ~/ 2;
      if (INCX == 1) {
        for (J = N; J >= 1; J--) {
          if (X[J] != Complex.zero) {
            if (NOUNIT) X[J] /= AP[KK];
            TEMP = X[J];
            K = KK - 1;
            for (I = J - 1; I >= 1; I--) {
              X[I] -= TEMP * AP[K];
              K--;
            }
          }
          KK -= J;
        }
      } else {
        JX = KX + (N - 1) * INCX;
        for (J = N; J >= 1; J--) {
          if (X[JX] != Complex.zero) {
            if (NOUNIT) X[JX] /= AP[KK];
            TEMP = X[JX];
            IX = JX;
            for (K = KK - 1; K >= KK - J + 1; K--) {
              IX -= INCX;
              X[IX] -= TEMP * AP[K];
            }
          }
          JX -= INCX;
          KK -= J;
        }
      }
    } else {
      KK = 1;
      if (INCX == 1) {
        for (J = 1; J <= N; J++) {
          if (X[J] != Complex.zero) {
            if (NOUNIT) X[J] /= AP[KK];
            TEMP = X[J];
            K = KK + 1;
            for (I = J + 1; I <= N; I++) {
              X[I] -= TEMP * AP[K];
              K++;
            }
          }
          KK += (N - J + 1);
        }
      } else {
        JX = KX;
        for (J = 1; J <= N; J++) {
          if (X[JX] != Complex.zero) {
            if (NOUNIT) X[JX] /= AP[KK];
            TEMP = X[JX];
            IX = JX;
            for (K = KK + 1; K <= KK + N - J; K++) {
              IX += INCX;
              X[IX] -= TEMP * AP[K];
            }
          }
          JX += INCX;
          KK += (N - J + 1);
        }
      }
    }
  } else {
    // Form  x := inv( A**T )*x  or  x := inv( A**H )*x.

    if (lsame(UPLO, 'U')) {
      KK = 1;
      if (INCX == 1) {
        for (J = 1; J <= N; J++) {
          TEMP = X[J];
          K = KK;
          if (NOCONJ) {
            for (I = 1; I <= J - 1; I++) {
              TEMP -= AP[K] * X[I];
              K++;
            }
            if (NOUNIT) TEMP /= AP[KK + J - 1];
          } else {
            for (I = 1; I <= J - 1; I++) {
              TEMP -= AP[K].conjugate() * X[I];
              K++;
            }
            if (NOUNIT) TEMP /= AP[KK + J - 1].conjugate();
          }
          X[J] = TEMP;
          KK += J;
        }
      } else {
        JX = KX;
        for (J = 1; J <= N; J++) {
          TEMP = X[JX];
          IX = KX;
          if (NOCONJ) {
            for (K = KK; K <= KK + J - 2; K++) {
              TEMP -= AP[K] * X[IX];
              IX += INCX;
            }
            if (NOUNIT) TEMP /= AP[KK + J - 1];
          } else {
            for (K = KK; K <= KK + J - 2; K++) {
              TEMP -= AP[K].conjugate() * X[IX];
              IX += INCX;
            }
            if (NOUNIT) TEMP /= AP[KK + J - 1].conjugate();
          }
          X[JX] = TEMP;
          JX += INCX;
          KK += J;
        }
      }
    } else {
      KK = (N * (N + 1)) ~/ 2;
      if (INCX == 1) {
        for (J = N; J >= 1; J--) {
          TEMP = X[J];
          K = KK;
          if (NOCONJ) {
            for (I = N; I >= J + 1; I--) {
              TEMP -= AP[K] * X[I];
              K--;
            }
            if (NOUNIT) TEMP /= AP[KK - N + J];
          } else {
            for (I = N; I >= J + 1; I--) {
              TEMP -= AP[K].conjugate() * X[I];
              K--;
            }
            if (NOUNIT) TEMP /= AP[KK - N + J].conjugate();
          }
          X[J] = TEMP;
          KK -= (N - J + 1);
        }
      } else {
        KX += (N - 1) * INCX;
        JX = KX;
        for (J = N; J >= 1; J--) {
          TEMP = X[JX];
          IX = KX;
          if (NOCONJ) {
            for (K = KK; K >= KK - (N - (J + 1)); K--) {
              TEMP -= AP[K] * X[IX];
              IX -= INCX;
            }
            if (NOUNIT) TEMP /= AP[KK - N + J];
          } else {
            for (K = KK; K >= KK - (N - (J + 1)); K--) {
              TEMP -= AP[K].conjugate() * X[IX];
              IX -= INCX;
            }
            if (NOUNIT) TEMP /= AP[KK - N + J].conjugate();
          }
          X[JX] = TEMP;
          JX -= INCX;
          KK -= (N - J + 1);
        }
      }
    }
  }
}
