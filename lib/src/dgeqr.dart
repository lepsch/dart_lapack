// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgeqrt.dart';
import 'package:dart_lapack/src/dlatsqr.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dgeqr(
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> T_,
  final int TSIZE,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final T = T_.having();
  final WORK = WORK_.having();
  bool LQUERY, LMINWS, MINT, MINW;
  int MB, NB, MINTSZ, NBLCKS, LWMIN, LWREQ;

  // Test the input arguments

  INFO.value = 0;

  LQUERY = (TSIZE == -1 || TSIZE == -2 || LWORK == -1 || LWORK == -2);

  MINT = false;
  MINW = false;
  if (TSIZE == -2 || LWORK == -2) {
    if (TSIZE != -1) MINT = true;
    if (LWORK != -1) MINW = true;
  }

  // Determine the block size

  if (min(M, N) > 0) {
    MB = ilaenv(1, 'DGEQR ', ' ', M, N, 1, -1);
    NB = ilaenv(1, 'DGEQR ', ' ', M, N, 2, -1);
  } else {
    MB = M;
    NB = 1;
  }
  if (MB > M || MB <= N) MB = M;
  if (NB > min(M, N) || NB < 1) NB = 1;
  MINTSZ = N + 5;
  if (MB > N && M > N) {
    if (((M - N) % (MB - N)) == 0) {
      NBLCKS = (M - N) ~/ (MB - N);
    } else {
      NBLCKS = (M - N) ~/ (MB - N) + 1;
    }
  } else {
    NBLCKS = 1;
  }

  // Determine if the workspace size satisfies minimal size

  LWMIN = max(1, N);
  LWREQ = max(1, N * NB);
  LMINWS = false;
  if ((TSIZE < max(1, NB * N * NBLCKS + 5) || LWORK < LWREQ) &&
      (LWORK >= N) &&
      (TSIZE >= MINTSZ) &&
      !LQUERY) {
    if (TSIZE < max(1, NB * N * NBLCKS + 5)) {
      LMINWS = true;
      NB = 1;
      MB = M;
    }
    if (LWORK < LWREQ) {
      LMINWS = true;
      NB = 1;
    }
  }

  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, M)) {
    INFO.value = -4;
  } else if (TSIZE < max(1, NB * N * NBLCKS + 5) && !LQUERY && !LMINWS) {
    INFO.value = -6;
  } else if ((LWORK < LWREQ) && !LQUERY && !LMINWS) {
    INFO.value = -8;
  }

  if (INFO.value == 0) {
    if (MINT) {
      T[1] = MINTSZ.toDouble();
    } else {
      T[1] = NB * N * NBLCKS + 5;
    }
    T[2] = MB.toDouble();
    T[3] = NB.toDouble();
    if (MINW) {
      WORK[1] = LWMIN.toDouble();
    } else {
      WORK[1] = LWREQ.toDouble();
    }
  }
  if (INFO.value != 0) {
    xerbla('DGEQR', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (min(M, N) == 0) {
    return;
  }

  // The QR Decomposition

  if ((M <= N) || (MB <= N) || (MB >= M)) {
    dgeqrt(M, N, NB, A, LDA, T(6).asMatrix(NB), NB, WORK, INFO);
  } else {
    dlatsqr(M, N, MB, NB, A, LDA, T(6).asMatrix(NB), NB, WORK, LWORK, INFO);
  }

  WORK[1] = LWREQ.toDouble();
}
