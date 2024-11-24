// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dspr.dart';
import 'package:lapack/src/blas/dtpsv.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dpptrf(
  final String UPLO,
  final int N,
  final Array<double> AP_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.having();
  const ONE = 1.0, ZERO = 0.0;
  bool UPPER;
  int J, JC, JJ;
  double AJJ;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  }
  if (INFO.value != 0) {
    xerbla('DPPTRF', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (UPPER) {
    // Compute the Cholesky factorization A = U**T*U.

    JJ = 0;
    for (J = 1; J <= N; J++) {
      JC = JJ + 1;
      JJ += J;

      // Compute elements 1:J-1 of column J.

      if (J > 1) dtpsv('Upper', 'Transpose', 'Non-unit', J - 1, AP, AP(JC), 1);

      // Compute U(J,J) and test for non-positive-definiteness.

      AJJ = AP[JJ] - ddot(J - 1, AP(JC), 1, AP(JC), 1);
      if (AJJ <= ZERO) {
        AP[JJ] = AJJ;
        INFO.value = J;
        return;
      }
      AP[JJ] = sqrt(AJJ);
    }
  } else {
    // Compute the Cholesky factorization A = L*L**T.

    JJ = 1;
    for (J = 1; J <= N; J++) {
      // Compute L(J,J) and test for non-positive-definiteness.

      AJJ = AP[JJ];
      if (AJJ <= ZERO) {
        AP[JJ] = AJJ;
        INFO.value = J;
        return;
      }
      AJJ = sqrt(AJJ);
      AP[JJ] = AJJ;

      // Compute elements J+1:N of column J and update the trailing
      // submatrix.

      if (J < N) {
        dscal(N - J, ONE / AJJ, AP(JJ + 1), 1);
        dspr('Lower', N - J, -ONE, AP(JJ + 1), 1, AP(JJ + N - J + 1));
        JJ += N - J + 1;
      }
    }
  }
}
