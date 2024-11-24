// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dtrsm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlaorhr_col_getrfnp2(
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> D_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final D = D_.having();
  const ONE = 1.0;
  double SFMIN;
  int I, N1, N2;
  final IINFO = Box(0);

  // Test the input parameters

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, M)) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('DLAORHR_COL_GETRFNP2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (min(M, N) == 0) return;

  if (M == 1) {
    // One row case, (also recursion termination case),
    // use unblocked code

    // Transfer the sign

    D[1] = -sign(ONE, A[1][1]);

    // Construct the row of U

    A[1][1] -= D[1];
  } else if (N == 1) {
    // One column case, (also recursion termination case),
    // use unblocked code

    // Transfer the sign

    D[1] = -sign(ONE, A[1][1]);

    // Construct the row of U

    A[1][1] -= D[1];

    // Scale the elements 2:M of the column

    // Determine machine safe minimum

    SFMIN = dlamch('S');

    // Construct the subdiagonal elements of L

    if (A[1][1].abs() >= SFMIN) {
      dscal(M - 1, ONE / A[1][1], A(2, 1).asArray(), 1);
    } else {
      for (I = 2; I <= M; I++) {
        A[I][1] /= A[1][1];
      }
    }
  } else {
    // Divide the matrix B into four submatrices

    N1 = min(M, N) ~/ 2;
    N2 = N - N1;

    // Factor B11, recursive call

    dlaorhr_col_getrfnp2(N1, N1, A, LDA, D, IINFO);

    // Solve for B21

    dtrsm('R', 'U', 'N', 'N', M - N1, N1, ONE, A, LDA, A(N1 + 1, 1), LDA);

    // Solve for B12

    dtrsm('L', 'L', 'N', 'U', N1, N2, ONE, A, LDA, A(1, N1 + 1), LDA);

    // Update B22, i.e. compute the Schur complement
    // B22 := B22 - B21*B12

    dgemm('N', 'N', M - N1, N2, N1, -ONE, A(N1 + 1, 1), LDA, A(1, N1 + 1), LDA,
        ONE, A(N1 + 1, N1 + 1), LDA);

    // Factor B22, recursive call

    dlaorhr_col_getrfnp2(M - N1, N2, A(N1 + 1, N1 + 1), LDA, D(N1 + 1), IINFO);
  }
}
