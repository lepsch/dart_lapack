// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zherk.dart';
import 'package:dart_lapack/src/blas/ztrsm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/disnan.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void zpotrf2(
  final String UPLO,
  final int N,
  final Matrix<Complex> A,
  final int LDA,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0, ZERO = 0.0;
  bool UPPER;
  int N1, N2;
  double AJJ;
  final IINFO = Box(0);

  // Test the input parameters

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('ZPOTRF2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // N=1 case

  if (N == 1) {
    // Test for non-positive-definiteness

    AJJ = A[1][1].real;
    if (AJJ <= ZERO || disnan(AJJ)) {
      INFO.value = 1;
      return;
    }

    // Factor

    A[1][1] = sqrt(AJJ).toComplex();

    // Use recursive code
  } else {
    N1 = N ~/ 2;
    N2 = N - N1;

    // Factor A11

    zpotrf2(UPLO, N1, A(1, 1), LDA, IINFO);
    if (IINFO.value != 0) {
      INFO.value = IINFO.value;
      return;
    }

    // Compute the Cholesky factorization A = U**H*U

    if (UPPER) {
      // Update and scale A12

      ztrsm('L', 'U', 'C', 'N', N1, N2, Complex.one, A(1, 1), LDA, A(1, N1 + 1),
          LDA);

      // Update and factor A22

      zherk(UPLO, 'C', N2, N1, -ONE, A(1, N1 + 1), LDA, ONE, A(N1 + 1, N1 + 1),
          LDA);
      zpotrf2(UPLO, N2, A(N1 + 1, N1 + 1), LDA, IINFO);
      if (IINFO.value != 0) {
        INFO.value = IINFO.value + N1;
        return;
      }

      // Compute the Cholesky factorization A = L*L**H
    } else {
      // Update and scale A21

      ztrsm('R', 'L', 'C', 'N', N2, N1, Complex.one, A(1, 1), LDA, A(N1 + 1, 1),
          LDA);

      // Update and factor A22

      zherk(UPLO, 'N', N2, N1, -ONE, A(N1 + 1, 1), LDA, ONE, A(N1 + 1, N1 + 1),
          LDA);
      zpotrf2(UPLO, N2, A(N1 + 1, N1 + 1), LDA, IINFO);
      if (IINFO.value != 0) {
        INFO.value = IINFO.value + N1;
        return;
      }
    }
  }
}
