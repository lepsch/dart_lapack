// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/dger.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgetf2(
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  const ONE = 1.0, ZERO = 0.0;
  double SFMIN;
  int I, J, JP;

  // Test the input parameters.

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, M)) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('DGETF2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) return;

  // Compute machine safe minimum

  SFMIN = dlamch('S');

  for (J = 1; J <= min(M, N); J++) {
    // Find pivot and test for singularity.

    JP = J - 1 + idamax(M - J + 1, A(J, J).asArray(), 1);
    IPIV[J] = JP;
    if (A[JP][J] != ZERO) {
      // Apply the interchange to columns 1:N.

      if (JP != J) dswap(N, A(J, 1).asArray(), LDA, A(JP, 1).asArray(), LDA);

      // Compute elements J+1:M of J-th column.

      if (J < M) {
        if (A[J][J].abs() >= SFMIN) {
          dscal(M - J, ONE / A[J][J], A(J + 1, J).asArray(), 1);
        } else {
          for (I = 1; I <= M - J; I++) {
            A[J + I][J] /= A[J][J];
          }
        }
      }
    } else if (INFO.value == 0) {
      INFO.value = J;
    }

    if (J < min(M, N)) {
      // Update trailing submatrix.

      dger(M - J, N - J, -ONE, A(J + 1, J).asArray(), 1, A(J, J + 1).asArray(),
          LDA, A(J + 1, J + 1), LDA);
    }
  }
}
