// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dtpqrt2.dart';
import 'package:dart_lapack/src/dtprfb.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dtpqrt(
  final int M,
  final int N,
  final int L,
  final int NB,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> T_,
  final int LDT,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final T = T_.having(ld: LDT);
  final WORK = WORK_.having();
  int I, IB, LB, MB;
  final IINFO = Box(0);

  // Test the input arguments

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (L < 0 || (L > min(M, N) && min(M, N) >= 0)) {
    INFO.value = -3;
  } else if (NB < 1 || (NB > N && N > 0)) {
    INFO.value = -4;
  } else if (LDA < max(1, N)) {
    INFO.value = -6;
  } else if (LDB < max(1, M)) {
    INFO.value = -8;
  } else if (LDT < NB) {
    INFO.value = -10;
  }
  if (INFO.value != 0) {
    xerbla('DTPQRT', -INFO.value);
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) return;

  for (I = 1; I <= N; I += NB) {
    // Compute the QR factorization of the current block

    IB = min(N - I + 1, NB);
    MB = min(M - L + I + IB - 1, M);
    if (I >= L) {
      LB = 0;
    } else {
      LB = MB - M + L - I + 1;
    }

    dtpqrt2(MB, IB, LB, A(I, I), LDA, B(1, I), LDB, T(1, I), LDT, IINFO);

    // Update by applying H**T to B(:,I+IB:N) from the left

    if (I + IB <= N) {
      dtprfb(
          'L',
          'T',
          'F',
          'C',
          MB,
          N - I - IB + 1,
          IB,
          LB,
          B(1, I),
          LDB,
          T(1, I),
          LDT,
          A(I, I + IB),
          LDA,
          B(1, I + IB),
          LDB,
          WORK.asMatrix(IB),
          IB);
    }
  }
}
