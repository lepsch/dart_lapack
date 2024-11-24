// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/blas/ztrsm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlaswp.dart';

void zgetrf2(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  double SFMIN;
  Complex TEMP;
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
    xerbla('ZGETRF2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) return;

  if (M == 1) {
    // Use unblocked code for one row case
    // Just need to handle IPIV and INFO

    IPIV[1] = 1;
    if (A[1][1] == Complex.zero) INFO.value = 1;
  } else if (N == 1) {
    // Use unblocked code for one column case

    // Compute machine safe minimum

    SFMIN = dlamch('S');

    // Find pivot and test for singularity

    I = izamax(M, A(1, 1).asArray(), 1);
    IPIV[1] = I;
    if (A[I][1] != Complex.zero) {
      // Apply the interchange

      if (I != 1) {
        TEMP = A[1][1];
        A[1][1] = A[I][1];
        A[I][1] = TEMP;
      }

      // Compute elements 2:M of the column

      if (A[1][1].abs() >= SFMIN) {
        zscal(M - 1, Complex.one / A[1][1], A(2, 1).asArray(), 1);
      } else {
        for (I = 1; I <= M - 1; I++) {
          A[1 + I][1] /= A[1][1];
        }
      }
    } else {
      INFO.value = 1;
    }
  } else {
    // Use recursive code

    N1 = min(M, N) ~/ 2;
    N2 = N - N1;

    //        [ A11 ]
    // Factor [ --- ]
    //        [ A21 ]

    zgetrf2(M, N1, A, LDA, IPIV, IINFO);
    if (INFO.value == 0 && IINFO.value > 0) INFO.value = IINFO.value;

    //                       [ A12 ]
    // Apply interchanges to [ --- ]
    //                       [ A22 ]

    zlaswp(N2, A(1, N1 + 1), LDA, 1, N1, IPIV, 1);

    // Solve A12

    ztrsm('L', 'L', 'N', 'U', N1, N2, Complex.one, A, LDA, A(1, N1 + 1), LDA);

    // Update A22

    zgemm('N', 'N', M - N1, N2, N1, -Complex.one, A(N1 + 1, 1), LDA,
        A(1, N1 + 1), LDA, Complex.one, A(N1 + 1, N1 + 1), LDA);

    // Factor A22

    zgetrf2(M - N1, N2, A(N1 + 1, N1 + 1), LDA, IPIV(N1 + 1), IINFO);

    // Adjust INFO and the pivot indices

    if (INFO.value == 0 && IINFO.value > 0) INFO.value = IINFO.value + N1;
    for (I = N1 + 1; I <= min(M, N); I++) {
      IPIV[I] += N1;
    }

    // Apply interchanges to A21

    zlaswp(N1, A(1, 1), LDA, N1 + 1, min(M, N), IPIV, 1);
  }
}
