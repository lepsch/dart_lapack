// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgeqrt.dart';
import 'package:dart_lapack/src/dtpqrt.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dlatsqr(
  final int M,
  final int N,
  final int MB,
  final int NB,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> T_,
  final int LDT,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final T = T_.having(ld: LDT);
  final WORK = WORK_.having();
  bool LQUERY;
  int I, II, KK, CTR, MINMN, LWMIN;

  // TEST THE INPUT ARGUMENTS

  INFO.value = 0;

  LQUERY = (LWORK == -1);

  MINMN = min(M, N);
  if (MINMN == 0) {
    LWMIN = 1;
  } else {
    LWMIN = N * NB;
  }

  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0 || M < N) {
    INFO.value = -2;
  } else if (MB < 1) {
    INFO.value = -3;
  } else if (NB < 1 || (NB > N && N > 0)) {
    INFO.value = -4;
  } else if (LDA < max(1, M)) {
    INFO.value = -6;
  } else if (LDT < NB) {
    INFO.value = -8;
  } else if (LWORK < LWMIN && !LQUERY) {
    INFO.value = -10;
  }

  if (INFO.value == 0) {
    WORK[1] = LWMIN.toDouble();
  }

  if (INFO.value != 0) {
    xerbla('DLATSQR', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (MINMN == 0) {
    return;
  }

  // The QR Decomposition

  if ((MB <= N) || (MB >= M)) {
    dgeqrt(M, N, NB, A, LDA, T, LDT, WORK, INFO);
    return;
  }

  KK = ((M - N) % (MB - N));
  II = M - KK + 1;

  // Compute the QR factorization of the first block A(1:MB,1:N)

  dgeqrt(MB, N, NB, A(1, 1), LDA, T, LDT, WORK, INFO);

  CTR = 1;
  for (I = MB + 1;
      (MB - N) < 0 ? I >= II - MB + N : I <= II - MB + N;
      I += (MB - N)) {
    // Compute the QR factorization of the current block A(I:I+MB-N,1:N)

    dtpqrt(MB - N, N, 0, NB, A(1, 1), LDA, A(I, 1), LDA, T(1, CTR * N + 1), LDT,
        WORK, INFO);
    CTR++;
  }

  // Compute the QR factorization of the last block A(II:M,1:N)

  if (II <= M) {
    dtpqrt(KK, N, 0, NB, A(1, 1), LDA, A(II, 1), LDA, T(1, CTR * N + 1), LDT,
        WORK, INFO);
  }

  WORK[1] = LWMIN.toDouble();
}
