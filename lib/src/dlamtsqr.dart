// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgemqrt.dart';
import 'package:dart_lapack/src/dtpmqrt.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dlamtsqr(
  final String SIDE,
  final String TRANS,
  final int M,
  final int N,
  final int K,
  final int MB,
  final int NB,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> T_,
  final int LDT,
  final Matrix<double> C_,
  final int LDC,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final T = T_.having(ld: LDT);
  final C = C_.having(ld: LDC);
  final WORK = WORK_.having();

  bool LEFT, RIGHT, TRAN, NOTRAN, LQUERY;
  int I, II, KK, LW, CTR, Q, MINMNK, LWMIN;

  // Test the input arguments

  INFO.value = 0;
  LQUERY = (LWORK == -1);
  NOTRAN = lsame(TRANS, 'N');
  TRAN = lsame(TRANS, 'T');
  LEFT = lsame(SIDE, 'L');
  RIGHT = lsame(SIDE, 'R');
  if (LEFT) {
    LW = N * NB;
    Q = M;
  } else {
    LW = MB * NB;
    Q = N;
  }

  MINMNK = min(M, min(N, K));
  if (MINMNK == 0) {
    LWMIN = 1;
  } else {
    LWMIN = max(1, LW);
  }

  if (!LEFT && !RIGHT) {
    INFO.value = -1;
  } else if (!TRAN && !NOTRAN) {
    INFO.value = -2;
  } else if (M < K) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (K < 0) {
    INFO.value = -5;
  } else if (K < NB || NB < 1) {
    INFO.value = -7;
  } else if (LDA < max(1, Q)) {
    INFO.value = -9;
  } else if (LDT < max(1, NB)) {
    INFO.value = -11;
  } else if (LDC < max(1, M)) {
    INFO.value = -13;
  } else if (LWORK < LWMIN && !LQUERY) {
    INFO.value = -15;
  }

  if (INFO.value == 0) {
    WORK[1] = LWMIN.toDouble();
  }

  if (INFO.value != 0) {
    xerbla('DLAMTSQR', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (MINMNK == 0) {
    return;
  }

  // Determine the block size if it is tall skinny or short and wide

  if ((MB <= K) || (MB >= max(M, max(N, K)))) {
    dgemqrt(SIDE, TRANS, M, N, K, NB, A, LDA, T, LDT, C, LDC, WORK, INFO);
    return;
  }

  if (LEFT && NOTRAN) {
    // Multiply Q to the last block of C

    KK = (M - K) % (MB - K);
    CTR = (M - K) ~/ (MB - K);
    if (KK > 0) {
      II = M - KK + 1;
      dtpmqrt('L', 'N', KK, N, K, 0, NB, A(II, 1), LDA, T(1, CTR * K + 1), LDT,
          C(1, 1), LDC, C(II, 1), LDC, WORK, INFO);
    } else {
      II = M + 1;
    }

    for (I = II - (MB - K); I >= MB + 1; I -= (MB - K)) {
      // Multiply Q to the current block of C (I:I+MB,1:N)

      CTR--;
      dtpmqrt('L', 'N', MB - K, N, K, 0, NB, A(I, 1), LDA, T(1, CTR * K + 1),
          LDT, C(1, 1), LDC, C(I, 1), LDC, WORK, INFO);
    }

    // Multiply Q to the first block of C (1:MB,1:N)

    dgemqrt(
        'L', 'N', MB, N, K, NB, A(1, 1), LDA, T, LDT, C(1, 1), LDC, WORK, INFO);
  } else if (LEFT && TRAN) {
    // Multiply Q to the first block of C

    KK = ((M - K) % (MB - K));
    II = M - KK + 1;
    CTR = 1;
    dgemqrt(
        'L', 'T', MB, N, K, NB, A(1, 1), LDA, T, LDT, C(1, 1), LDC, WORK, INFO);

    for (I = MB + 1;
        (MB - K) < 0 ? I >= II - MB + K : I <= II - MB + K;
        I += (MB - K)) {
      // Multiply Q to the current block of C (I:I+MB,1:N)

      dtpmqrt('L', 'T', MB - K, N, K, 0, NB, A(I, 1), LDA, T(1, CTR * K + 1),
          LDT, C(1, 1), LDC, C(I, 1), LDC, WORK, INFO);
      CTR++;
    }
    if (II <= M) {
      // Multiply Q to the last block of C

      dtpmqrt('L', 'T', KK, N, K, 0, NB, A(II, 1), LDA, T(1, CTR * K + 1), LDT,
          C(1, 1), LDC, C(II, 1), LDC, WORK, INFO);
    }
  } else if (RIGHT && TRAN) {
    // Multiply Q to the last block of C

    KK = (N - K) % (MB - K);
    CTR = (N - K) ~/ (MB - K);
    if (KK > 0) {
      II = N - KK + 1;
      dtpmqrt('R', 'T', M, KK, K, 0, NB, A(II, 1), LDA, T(1, CTR * K + 1), LDT,
          C(1, 1), LDC, C(1, II), LDC, WORK, INFO);
    } else {
      II = N + 1;
    }

    for (I = II - (MB - K); I >= MB + 1; I -= (MB - K)) {
      // Multiply Q to the current block of C (1:M,I:I+MB)

      CTR--;
      dtpmqrt('R', 'T', M, MB - K, K, 0, NB, A(I, 1), LDA, T(1, CTR * K + 1),
          LDT, C(1, 1), LDC, C(1, I), LDC, WORK, INFO);
    }

    // Multiply Q to the first block of C (1:M,1:MB)

    dgemqrt(
        'R', 'T', M, MB, K, NB, A(1, 1), LDA, T, LDT, C(1, 1), LDC, WORK, INFO);
  } else if (RIGHT && NOTRAN) {
    // Multiply Q to the first block of C

    KK = ((N - K) % (MB - K));
    II = N - KK + 1;
    CTR = 1;
    dgemqrt(
        'R', 'N', M, MB, K, NB, A(1, 1), LDA, T, LDT, C(1, 1), LDC, WORK, INFO);

    for (I = MB + 1;
        (MB - K) < 0 ? I >= II - MB + K : I <= II - MB + K;
        I += (MB - K)) {
      // Multiply Q to the current block of C (1:M,I:I+MB)

      dtpmqrt('R', 'N', M, MB - K, K, 0, NB, A(I, 1), LDA, T(1, CTR * K + 1),
          LDT, C(1, 1), LDC, C(1, I), LDC, WORK, INFO);
      CTR++;
    }
    if (II <= N) {
      // Multiply Q to the last block of C

      dtpmqrt('R', 'N', M, KK, K, 0, NB, A(II, 1), LDA, T(1, CTR * K + 1), LDT,
          C(1, 1), LDC, C(1, II), LDC, WORK, INFO);
    }
  }

  WORK[1] = LWMIN.toDouble();
}
