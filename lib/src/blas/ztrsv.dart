// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void ztrsv(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> X_,
  final int INCX,
) {
// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final X = X_.having();
  Complex TEMP;
  int I, INFO, IX, J, JX, KX = 0;
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
  } else if (LDA < max(1, N)) {
    INFO = 6;
  } else if (INCX == 0) {
    INFO = 8;
  }
  if (INFO != 0) {
    xerbla('ZTRSV', INFO);
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

  // Start the operations. In this version the elements of A are
  // accessed sequentially with one pass through A.

  if (lsame(TRANS, 'N')) {
    // Form  x := inv( A )*x.

    if (lsame(UPLO, 'U')) {
      if (INCX == 1) {
        for (J = N; J >= 1; J--) {
          if (X[J] != Complex.zero) {
            if (NOUNIT) X[J] /= A[J][J];
            TEMP = X[J];
            for (I = J - 1; I >= 1; I--) {
              X[I] -= TEMP * A[I][J];
            }
          }
        }
      } else {
        JX = KX + (N - 1) * INCX;
        for (J = N; J >= 1; J--) {
          if (X[JX] != Complex.zero) {
            if (NOUNIT) X[JX] /= A[J][J];
            TEMP = X[JX];
            IX = JX;
            for (I = J - 1; I >= 1; I--) {
              IX -= INCX;
              X[IX] -= TEMP * A[I][J];
            }
          }
          JX -= INCX;
        }
      }
    } else {
      if (INCX == 1) {
        for (J = 1; J <= N; J++) {
          if (X[J] != Complex.zero) {
            if (NOUNIT) X[J] /= A[J][J];
            TEMP = X[J];
            for (I = J + 1; I <= N; I++) {
              X[I] -= TEMP * A[I][J];
            }
          }
        }
      } else {
        JX = KX;
        for (J = 1; J <= N; J++) {
          if (X[JX] != Complex.zero) {
            if (NOUNIT) X[JX] /= A[J][J];
            TEMP = X[JX];
            IX = JX;
            for (I = J + 1; I <= N; I++) {
              IX += INCX;
              X[IX] -= TEMP * A[I][J];
            }
          }
          JX += INCX;
        }
      }
    }
  } else {
    // Form  x := inv( A**T )*x  or  x := inv( A**H )*x.

    if (lsame(UPLO, 'U')) {
      if (INCX == 1) {
        for (J = 1; J <= N; J++) {
          TEMP = X[J];
          if (NOCONJ) {
            for (I = 1; I <= J - 1; I++) {
              TEMP -= A[I][J] * X[I];
            }
            if (NOUNIT) TEMP /= A[J][J];
          } else {
            for (I = 1; I <= J - 1; I++) {
              TEMP -= A[I][J].conjugate() * X[I];
            }
            if (NOUNIT) TEMP /= A[J][J].conjugate();
          }
          X[J] = TEMP;
        }
      } else {
        JX = KX;
        for (J = 1; J <= N; J++) {
          IX = KX;
          TEMP = X[JX];
          if (NOCONJ) {
            for (I = 1; I <= J - 1; I++) {
              TEMP -= A[I][J] * X[IX];
              IX += INCX;
            }
            if (NOUNIT) TEMP /= A[J][J];
          } else {
            for (I = 1; I <= J - 1; I++) {
              TEMP -= A[I][J].conjugate() * X[IX];
              IX += INCX;
            }
            if (NOUNIT) TEMP /= A[J][J].conjugate();
          }
          X[JX] = TEMP;
          JX += INCX;
        }
      }
    } else {
      if (INCX == 1) {
        for (J = N; J >= 1; J--) {
          TEMP = X[J];
          if (NOCONJ) {
            for (I = N; I >= J + 1; I--) {
              TEMP -= A[I][J] * X[I];
            }
            if (NOUNIT) TEMP /= A[J][J];
          } else {
            for (I = N; I >= J + 1; I--) {
              TEMP -= A[I][J].conjugate() * X[I];
            }
            if (NOUNIT) TEMP /= A[J][J].conjugate();
          }
          X[J] = TEMP;
        }
      } else {
        KX += (N - 1) * INCX;
        JX = KX;
        for (J = N; J >= 1; J--) {
          IX = KX;
          TEMP = X[JX];
          if (NOCONJ) {
            for (I = N; I >= J + 1; I--) {
              TEMP -= A[I][J] * X[IX];
              IX -= INCX;
            }
            if (NOUNIT) TEMP /= A[J][J];
          } else {
            for (I = N; I >= J + 1; I--) {
              TEMP -= A[I][J].conjugate() * X[IX];
              IX -= INCX;
            }
            if (NOUNIT) TEMP /= A[J][J].conjugate();
          }
          X[JX] = TEMP;
          JX -= INCX;
        }
      }
    }
  }
}
