// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zgeqrt2.dart';
import 'package:dart_lapack/src/zgeqrt3.dart';
import 'package:dart_lapack/src/zlarfb.dart';

void zgeqrt(
  final int M,
  final int N,
  final int NB,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> T_,
  final int LDT,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final T = T_.having(ld: LDT);
  final WORK = WORK_.having();
  int I, IB, K;
  const USE_RECURSIVE_QR = true;
  final IINFO = Box(0);

  // Test the input arguments

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NB < 1 || (NB > min(M, N) && min(M, N) > 0)) {
    INFO.value = -3;
  } else if (LDA < max(1, M)) {
    INFO.value = -5;
  } else if (LDT < NB) {
    INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('ZGEQRT', -INFO.value);
    return;
  }

  // Quick return if possible

  K = min(M, N);
  if (K == 0) return;

  // Blocked loop of length K

  for (I = 1; I <= K; I += NB) {
    IB = min(K - I + 1, NB);

    // Compute the QR factorization of the current block A(I:M,I:I+IB-1)

    if (USE_RECURSIVE_QR) {
      zgeqrt3(M - I + 1, IB, A(I, I), LDA, T(1, I), LDT, IINFO);
      // ignore: dead_code
    } else {
      zgeqrt2(M - I + 1, IB, A(I, I), LDA, T(1, I), LDT, IINFO);
    }
    if (I + IB <= N) {
      // Update by applying H**H to A(I:M,I+IB:N) from the left

      zlarfb(
          'L',
          'C',
          'F',
          'C',
          M - I + 1,
          N - I - IB + 1,
          IB,
          A(I, I),
          LDA,
          T(1, I),
          LDT,
          A(I, I + IB),
          LDA,
          WORK.asMatrix(N - I - IB + 1),
          N - I - IB + 1);
    }
  }
}
