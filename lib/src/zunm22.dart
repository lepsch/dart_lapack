// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/blas/ztrmm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlacpy.dart';

void zunm22(
  final String SIDE,
  final String TRANS,
  final int M,
  final int N,
  final int N1,
  final int N2,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Matrix<Complex> C_,
  final int LDC,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final Q = Q_.having(ld: LDQ);
  final C = C_.having(ld: LDC);
  final WORK = WORK_.having();
  bool LEFT, LQUERY, NOTRAN;
  int I, LDWORK, LEN, LWKOPT = 0, NB, NQ, NW;

  // Test the input arguments

  INFO.value = 0;
  LEFT = lsame(SIDE, 'L');
  NOTRAN = lsame(TRANS, 'N');
  LQUERY = (LWORK == -1);

  // NQ is the order of Q;
  // NW is the minimum dimension of WORK.

  if (LEFT) {
    NQ = M;
  } else {
    NQ = N;
  }
  NW = NQ;
  if (N1 == 0 || N2 == 0) NW = 1;
  if (!LEFT && !lsame(SIDE, 'R')) {
    INFO.value = -1;
  } else if (!lsame(TRANS, 'N') && !lsame(TRANS, 'C')) {
    INFO.value = -2;
  } else if (M < 0) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (N1 < 0 || N1 + N2 != NQ) {
    INFO.value = -5;
  } else if (N2 < 0) {
    INFO.value = -6;
  } else if (LDQ < max(1, NQ)) {
    INFO.value = -8;
  } else if (LDC < max(1, M)) {
    INFO.value = -10;
  } else if (LWORK < NW && !LQUERY) {
    INFO.value = -12;
  }

  if (INFO.value == 0) {
    LWKOPT = M * N;
    WORK[1] = LWKOPT.toComplex();
  }

  if (INFO.value != 0) {
    xerbla('ZUNM22', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) {
    WORK[1] = Complex.one;
    return;
  }

  // Degenerate cases (N1 = 0 or N2 = 0) are handled using ZTRMM.

  if (N1 == 0) {
    ztrmm(SIDE, 'Upper', TRANS, 'Non-Unit', M, N, Complex.one, Q, LDQ, C, LDC);
    WORK[1] = Complex.one;
    return;
  } else if (N2 == 0) {
    ztrmm(SIDE, 'Lower', TRANS, 'Non-Unit', M, N, Complex.one, Q, LDQ, C, LDC);
    WORK[1] = Complex.one;
    return;
  }

  // Compute the largest chunk size available from the workspace.

  NB = max(1, min(LWORK, LWKOPT) ~/ NQ);

  if (LEFT) {
    if (NOTRAN) {
      for (I = 1; I <= N; I += NB) {
        LEN = min(NB, N - I + 1);
        LDWORK = M;

        // Multiply bottom part of C by Q12.

        zlacpy(
            'All', N1, LEN, C(N2 + 1, I), LDC, WORK.asMatrix(LDWORK), LDWORK);
        ztrmm('Left', 'Lower', 'No Transpose', 'Non-Unit', N1, LEN, Complex.one,
            Q(1, N2 + 1), LDQ, WORK.asMatrix(LDWORK), LDWORK);

        // Multiply top part of C by Q11.

        zgemm('No Transpose', 'No Transpose', N1, LEN, N2, Complex.one, Q, LDQ,
            C(1, I), LDC, Complex.one, WORK.asMatrix(LDWORK), LDWORK);

        // Multiply top part of C by Q21.

        zlacpy('All', N2, LEN, C(1, I), LDC, WORK(N1 + 1).asMatrix(LDWORK),
            LDWORK);
        ztrmm('Left', 'Upper', 'No Transpose', 'Non-Unit', N2, LEN, Complex.one,
            Q(N1 + 1, 1), LDQ, WORK(N1 + 1).asMatrix(LDWORK), LDWORK);

        // Multiply bottom part of C by Q22.

        zgemm(
            'No Transpose',
            'No Transpose',
            N2,
            LEN,
            N1,
            Complex.one,
            Q(N1 + 1, N2 + 1),
            LDQ,
            C(N2 + 1, I),
            LDC,
            Complex.one,
            WORK(N1 + 1).asMatrix(LDWORK),
            LDWORK);

        // Copy everything back.

        zlacpy('All', M, LEN, WORK.asMatrix(LDWORK), LDWORK, C(1, I), LDC);
      }
    } else {
      for (I = 1; I <= N; I += NB) {
        LEN = min(NB, N - I + 1);
        LDWORK = M;

        // Multiply bottom part of C by Q21**H.

        zlacpy(
            'All', N2, LEN, C(N1 + 1, I), LDC, WORK.asMatrix(LDWORK), LDWORK);
        ztrmm('Left', 'Upper', 'Conjugate', 'Non-Unit', N2, LEN, Complex.one,
            Q(N1 + 1, 1), LDQ, WORK.asMatrix(LDWORK), LDWORK);

        // Multiply top part of C by Q11**H.

        zgemm('Conjugate', 'No Transpose', N2, LEN, N1, Complex.one, Q, LDQ,
            C(1, I), LDC, Complex.one, WORK.asMatrix(LDWORK), LDWORK);

        // Multiply top part of C by Q12**H.

        zlacpy('All', N1, LEN, C(1, I), LDC, WORK(N2 + 1).asMatrix(LDWORK),
            LDWORK);
        ztrmm('Left', 'Lower', 'Conjugate', 'Non-Unit', N1, LEN, Complex.one,
            Q(1, N2 + 1), LDQ, WORK(N2 + 1).asMatrix(LDWORK), LDWORK);

        // Multiply bottom part of C by Q22**H.

        zgemm(
            'Conjugate',
            'No Transpose',
            N1,
            LEN,
            N2,
            Complex.one,
            Q(N1 + 1, N2 + 1),
            LDQ,
            C(N1 + 1, I),
            LDC,
            Complex.one,
            WORK(N2 + 1).asMatrix(LDWORK),
            LDWORK);

        // Copy everything back.

        zlacpy('All', M, LEN, WORK.asMatrix(LDWORK), LDWORK, C(1, I), LDC);
      }
    }
  } else {
    if (NOTRAN) {
      for (I = 1; I <= M; I += NB) {
        LEN = min(NB, M - I + 1);
        LDWORK = LEN;

        // Multiply right part of C by Q21.

        zlacpy(
            'All', LEN, N2, C(I, N1 + 1), LDC, WORK.asMatrix(LDWORK), LDWORK);
        ztrmm('Right', 'Upper', 'No Transpose', 'Non-Unit', LEN, N2,
            Complex.one, Q(N1 + 1, 1), LDQ, WORK.asMatrix(LDWORK), LDWORK);

        // Multiply left part of C by Q11.

        zgemm('No Transpose', 'No Transpose', LEN, N2, N1, Complex.one, C(I, 1),
            LDC, Q, LDQ, Complex.one, WORK.asMatrix(LDWORK), LDWORK);

        // Multiply left part of C by Q12.

        zlacpy('All', LEN, N1, C(I, 1), LDC,
            WORK(1 + N2 * LDWORK).asMatrix(LDWORK), LDWORK);
        ztrmm(
            'Right',
            'Lower',
            'No Transpose',
            'Non-Unit',
            LEN,
            N1,
            Complex.one,
            Q(1, N2 + 1),
            LDQ,
            WORK(1 + N2 * LDWORK).asMatrix(LDWORK),
            LDWORK);

        // Multiply right part of C by Q22.

        zgemm(
            'No Transpose',
            'No Transpose',
            LEN,
            N1,
            N2,
            Complex.one,
            C(I, N1 + 1),
            LDC,
            Q(N1 + 1, N2 + 1),
            LDQ,
            Complex.one,
            WORK(1 + N2 * LDWORK).asMatrix(LDWORK),
            LDWORK);

        // Copy everything back.

        zlacpy('All', LEN, N, WORK.asMatrix(LDWORK), LDWORK, C(I, 1), LDC);
      }
    } else {
      for (I = 1; I <= M; I += NB) {
        LEN = min(NB, M - I + 1);
        LDWORK = LEN;

        // Multiply right part of C by Q12**H.

        zlacpy(
            'All', LEN, N1, C(I, N2 + 1), LDC, WORK.asMatrix(LDWORK), LDWORK);
        ztrmm('Right', 'Lower', 'Conjugate', 'Non-Unit', LEN, N1, Complex.one,
            Q(1, N2 + 1), LDQ, WORK.asMatrix(LDWORK), LDWORK);

        // Multiply left part of C by Q11**H.

        zgemm('No Transpose', 'Conjugate', LEN, N1, N2, Complex.one, C(I, 1),
            LDC, Q, LDQ, Complex.one, WORK.asMatrix(LDWORK), LDWORK);

        // Multiply left part of C by Q21**H.

        zlacpy('All', LEN, N2, C(I, 1), LDC,
            WORK(1 + N1 * LDWORK).asMatrix(LDWORK), LDWORK);
        ztrmm('Right', 'Upper', 'Conjugate', 'Non-Unit', LEN, N2, Complex.one,
            Q(N1 + 1, 1), LDQ, WORK(1 + N1 * LDWORK).asMatrix(LDWORK), LDWORK);

        // Multiply right part of C by Q22**H.

        zgemm(
            'No Transpose',
            'Conjugate',
            LEN,
            N2,
            N1,
            Complex.one,
            C(I, N2 + 1),
            LDC,
            Q(N1 + 1, N2 + 1),
            LDQ,
            Complex.one,
            WORK(1 + N1 * LDWORK).asMatrix(LDWORK),
            LDWORK);

        // Copy everything back.

        zlacpy('All', LEN, N, WORK.asMatrix(LDWORK), LDWORK, C(I, 1), LDC);
      }
    }
  }

  WORK[1] = LWKOPT.toComplex();
}
