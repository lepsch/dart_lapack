// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zher2(
  final String UPLO,
  final int N,
  final Complex ALPHA,
  final Array<Complex> X_,
  final int INCX,
  final Array<Complex> Y_,
  final int INCY,
  final Matrix<Complex> A_,
  final int LDA,
) {
// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having();
  final Y = Y_.having();
  final A = A_.having(ld: LDA);
  Complex TEMP1, TEMP2;
  int I, INFO, IX, IY, J, JX = 0, JY = 0, KX = 0, KY = 0;

  // Test the input parameters.

  INFO = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO = 1;
  } else if (N < 0) {
    INFO = 2;
  } else if (INCX == 0) {
    INFO = 5;
  } else if (INCY == 0) {
    INFO = 7;
  } else if (LDA < max(1, N)) {
    INFO = 9;
  }
  if (INFO != 0) {
    xerbla('ZHER2', INFO);
    return;
  }

  // Quick return if possible.

  if ((N == 0) || (ALPHA == Complex.zero)) return;

  // Set up the start points in X and Y if the increments are not both
  // unity.

  if ((INCX != 1) || (INCY != 1)) {
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
    JX = KX;
    JY = KY;
  }

  // Start the operations. In this version the elements of A are
  // accessed sequentially with one pass through the triangular part
  // of A.

  if (lsame(UPLO, 'U')) {
    // Form  A  when A is stored in the upper triangle.

    if ((INCX == 1) && (INCY == 1)) {
      for (J = 1; J <= N; J++) {
        if ((X[J] != Complex.zero) || (Y[J] != Complex.zero)) {
          TEMP1 = ALPHA * Y[J].conjugate();
          TEMP2 = (ALPHA * X[J]).conjugate();
          for (I = 1; I <= J - 1; I++) {
            A[I][J] += X[I] * TEMP1 + Y[I] * TEMP2;
          }
          A[J][J] = Complex(A[J][J].real + (X[J] * TEMP1 + Y[J] * TEMP2).real);
        } else {
          A[J][J] = A[J][J].real.toComplex();
        }
      }
    } else {
      for (J = 1; J <= N; J++) {
        if ((X[JX] != Complex.zero) || (Y[JY] != Complex.zero)) {
          TEMP1 = ALPHA * Y[JY].conjugate();
          TEMP2 = (ALPHA * X[JX]).conjugate();
          IX = KX;
          IY = KY;
          for (I = 1; I <= J - 1; I++) {
            A[I][J] += X[IX] * TEMP1 + Y[IY] * TEMP2;
            IX += INCX;
            IY += INCY;
          }
          A[J][J] =
              (A[J][J].real + (X[JX] * TEMP1 + Y[JY] * TEMP2).real).toComplex();
        } else {
          A[J][J] = A[J][J].real.toComplex();
        }
        JX += INCX;
        JY += INCY;
      }
    }
  } else {
    // Form  A  when A is stored in the lower triangle.

    if ((INCX == 1) && (INCY == 1)) {
      for (J = 1; J <= N; J++) {
        if ((X[J] != Complex.zero) || (Y[J] != Complex.zero)) {
          TEMP1 = ALPHA * Y[J].conjugate();
          TEMP2 = (ALPHA * X[J]).conjugate();
          A[J][J] =
              (A[J][J].real + (X[J] * TEMP1 + Y[J] * TEMP2).real).toComplex();
          for (I = J + 1; I <= N; I++) {
            A[I][J] += X[I] * TEMP1 + Y[I] * TEMP2;
          }
        } else {
          A[J][J] = A[J][J].real.toComplex();
        }
      }
    } else {
      for (J = 1; J <= N; J++) {
        if ((X[JX] != Complex.zero) || (Y[JY] != Complex.zero)) {
          TEMP1 = ALPHA * Y[JY].conjugate();
          TEMP2 = (ALPHA * X[JX]).conjugate();
          A[J][J] =
              (A[J][J].real + (X[JX] * TEMP1 + Y[JY] * TEMP2).real).toComplex();
          IX = JX;
          IY = JY;
          for (I = J + 1; I <= N; I++) {
            IX += INCX;
            IY += INCY;
            A[I][J] += X[IX] * TEMP1 + Y[IY] * TEMP2;
          }
        } else {
          A[J][J] = A[J][J].real.toComplex();
        }
        JX += INCX;
        JY += INCY;
      }
    }
  }
}
