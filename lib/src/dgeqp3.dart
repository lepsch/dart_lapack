// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dnrm2.dart';
import 'package:dart_lapack/src/blas/dswap.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgeqrf.dart';
import 'package:dart_lapack/src/dlaqp2.dart';
import 'package:dart_lapack/src/dlaqps.dart';
import 'package:dart_lapack/src/dormqr.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dgeqp3(
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> JPVT_,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final JPVT = JPVT_.having();
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  const INB = 1, INBMIN = 2, IXOVER = 3;
  bool LQUERY;
  int IWS = 0,
      J,
      JB,
      LWKOPT,
      MINMN = 0,
      MINWS,
      NA,
      NB,
      NBMIN,
      NFXD,
      NX,
      SM,
      SMINMN,
      SN,
      TOPBMN = 0;
  final FJB = Box(0);

  // Test input arguments
  INFO.value = 0;
  LQUERY = (LWORK == -1);
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, M)) {
    INFO.value = -4;
  }

  if (INFO.value == 0) {
    MINMN = min(M, N);
    if (MINMN == 0) {
      IWS = 1;
      LWKOPT = 1;
    } else {
      IWS = 3 * N + 1;
      NB = ilaenv(INB, 'DGEQRF', ' ', M, N, -1, -1);
      LWKOPT = 2 * N + (N + 1) * NB;
    }
    WORK[1] = LWKOPT.toDouble();

    if ((LWORK < IWS) && !LQUERY) {
      INFO.value = -8;
    }
  }

  if (INFO.value != 0) {
    xerbla('DGEQP3', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Move initial columns up front.
  NFXD = 1;
  for (J = 1; J <= N; J++) {
    if (JPVT[J] != 0) {
      if (J != NFXD) {
        dswap(M, A(1, J).asArray(), 1, A(1, NFXD).asArray(), 1);
        JPVT[J] = JPVT[NFXD];
        JPVT[NFXD] = J;
      } else {
        JPVT[J] = J;
      }
      NFXD++;
    } else {
      JPVT[J] = J;
    }
  }
  NFXD--;

  // Factorize fixed columns

  // Compute the QR factorization of fixed columns and update
  // remaining columns.
  if (NFXD > 0) {
    NA = min(M, NFXD);
    // CALL DGEQR2( M, NA, A, LDA, TAU, WORK, INFO )
    dgeqrf(M, NA, A, LDA, TAU, WORK, LWORK, INFO);
    IWS = max(IWS, WORK[1].toInt());
    if (NA < N) {
      // CALL dorm2r( 'Left', 'Transpose', M, N-NA, NA, A, LDA, TAU, A( 1, NA+1 ), LDA, WORK, INFO )
      dormqr('Left', 'Transpose', M, N - NA, NA, A, LDA, TAU, A(1, NA + 1), LDA,
          WORK, LWORK, INFO);
      IWS = max(IWS, WORK[1].toInt());
    }
  }

  // Factorize free columns
  if (NFXD < MINMN) {
    SM = M - NFXD;
    SN = N - NFXD;
    SMINMN = MINMN - NFXD;

    // Determine the block size.
    NB = ilaenv(INB, 'DGEQRF', ' ', SM, SN, -1, -1);
    NBMIN = 2;
    NX = 0;

    if ((NB > 1) && (NB < SMINMN)) {
      // Determine when to cross over from blocked to unblocked code.
      NX = max(0, ilaenv(IXOVER, 'DGEQRF', ' ', SM, SN, -1, -1));

      if (NX < SMINMN) {
        // Determine if workspace is large enough for blocked code.
        MINWS = 2 * SN + (SN + 1) * NB;
        IWS = max(IWS, MINWS);
        if (LWORK < MINWS) {
          // Not enough workspace to use optimal NB: Reduce NB and
          // determine the minimum value of NB.
          NB = (LWORK - 2 * SN) ~/ (SN + 1);
          NBMIN = max(2, ilaenv(INBMIN, 'DGEQRF', ' ', SM, SN, -1, -1));
        }
      }
    }

    // Initialize partial column norms. The first N elements of work
    // store the exact column norms.
    for (J = NFXD + 1; J <= N; J++) {
      WORK[J] = dnrm2(SM, A(NFXD + 1, J).asArray(), 1);
      WORK[N + J] = WORK[J];
    }

    if ((NB >= NBMIN) && (NB < SMINMN) && (NX < SMINMN)) {
      // Use blocked code initially.
      J = NFXD + 1;

      // Compute factorization: while loop.
      TOPBMN = MINMN - NX;
      while (true) {
        if (J <= TOPBMN) {
          JB = min(NB, TOPBMN - J + 1);

          // Factorize JB columns among columns J:N.
          dlaqps(
            M,
            N - J + 1,
            J - 1,
            JB,
            FJB,
            A(1, J),
            LDA,
            JPVT(J),
            TAU(J),
            WORK(J),
            WORK(N + J),
            WORK(2 * N + 1),
            WORK(2 * N + JB + 1).asMatrix(N - J + 1),
            N - J + 1,
          );

          J += FJB.value;
          continue;
        }
        break;
      }
    } else {
      J = NFXD + 1;
    }

    // Use unblocked code to factor the last or only block.
    if (J <= MINMN) {
      dlaqp2(M, N - J + 1, J - 1, A(1, J), LDA, JPVT(J), TAU(J), WORK(J),
          WORK(N + J), WORK(2 * N + 1));
    }
  }

  WORK[1] = IWS.toDouble();
}
