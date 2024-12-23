// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dsyrk.dart';
import 'package:dart_lapack/src/blas/dtrsm.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/disnan.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dpotrf2(
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  const ONE = 1.0, ZERO = 0.0;
  bool UPPER;
  int N1, N2;
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
    xerbla('DPOTRF2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // N=1 case

  if (N == 1) {
    // Test for non-positive-definiteness

    if (A[1][1] <= ZERO || disnan(A[1][1])) {
      INFO.value = 1;
      return;
    }

    // Factor

    A[1][1] = sqrt(A[1][1]);

    // Use recursive code
  } else {
    N1 = N ~/ 2;
    N2 = N - N1;

    // Factor A11

    dpotrf2(UPLO, N1, A(1, 1), LDA, IINFO);
    if (IINFO.value != 0) {
      INFO.value = IINFO.value;
      return;
    }

    // Compute the Cholesky factorization A = U**T*U

    if (UPPER) {
      // Update and scale A12

      dtrsm('L', 'U', 'T', 'N', N1, N2, ONE, A(1, 1), LDA, A(1, N1 + 1), LDA);

      // Update and factor A22

      dsyrk(UPLO, 'T', N2, N1, -ONE, A(1, N1 + 1), LDA, ONE, A(N1 + 1, N1 + 1),
          LDA);
      dpotrf2(UPLO, N2, A(N1 + 1, N1 + 1), LDA, IINFO);
      if (IINFO.value != 0) {
        INFO.value = IINFO.value + N1;
        return;
      }

      // Compute the Cholesky factorization A = L*L**T
    } else {
      // Update and scale A21

      dtrsm('R', 'L', 'T', 'N', N2, N1, ONE, A(1, 1), LDA, A(N1 + 1, 1), LDA);

      // Update and factor A22

      dsyrk(UPLO, 'N', N2, N1, -ONE, A(N1 + 1, 1), LDA, ONE, A(N1 + 1, N1 + 1),
          LDA);
      dpotrf2(UPLO, N2, A(N1 + 1, N1 + 1), LDA, IINFO);
      if (IINFO.value != 0) {
        INFO.value = IINFO.value + N1;
        return;
      }
    }
  }
}
