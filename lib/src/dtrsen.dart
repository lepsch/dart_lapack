// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dtrexc.dart';
import 'package:dart_lapack/src/dtrsyl.dart';
import 'package:dart_lapack/src/dlacn2.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dtrsen(
  final String JOB,
  final String COMPQ,
  final Array<bool> SELECT_,
  final int N,
  final Matrix<double> T_,
  final int LDT,
  final Matrix<double> Q_,
  final int LDQ,
  final Array<double> WR_,
  final Array<double> WI_,
  final Box<int> M,
  final Box<double> S,
  final Box<double> SEP,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Box<int> INFO,
) {
  final SELECT = SELECT_.having();
  final T = T_.having(ld: LDT);
  final Q = Q_.having(ld: LDQ);
  final WR = WR_.having();
  final WI = WI_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool LQUERY, PAIR, SWAP, WANTBH, WANTQ, WANTS, WANTSP;
  int K, LIWMIN = 0, LWMIN = 0, N1 = 0, N2 = 0, NN = 0;
  double RNORM = 0;
  final ISAVE = Array<int>(3);
  final IERR = Box(0), KK = Box(0), KS = Box(0), KASE = Box(0);
  final SCALE = Box(0.0), EST = Box(0.0);

  // Decode and test the input parameters

  WANTBH = lsame(JOB, 'B');
  WANTS = lsame(JOB, 'E') || WANTBH;
  WANTSP = lsame(JOB, 'V') || WANTBH;
  WANTQ = lsame(COMPQ, 'V');

  INFO.value = 0;
  LQUERY = (LWORK == -1);
  if (!lsame(JOB, 'N') && !WANTS && !WANTSP) {
    INFO.value = -1;
  } else if (!lsame(COMPQ, 'N') && !WANTQ) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (LDT < max(1, N)) {
    INFO.value = -6;
  } else if (LDQ < 1 || (WANTQ && LDQ < N)) {
    INFO.value = -8;
  } else {
    // Set M to the dimension of the specified invariant subspace,
    // and test LWORK and LIWORK.

    M.value = 0;
    PAIR = false;
    for (K = 1; K <= N; K++) {
      if (PAIR) {
        PAIR = false;
      } else {
        if (K < N) {
          if (T[K + 1][K] == ZERO) {
            if (SELECT[K]) M.value++;
          } else {
            PAIR = true;
            if (SELECT[K] || SELECT[K + 1]) M.value += 2;
          }
        } else {
          if (SELECT[N]) M.value++;
        }
      }
    }

    N1 = M.value;
    N2 = N - M.value;
    NN = N1 * N2;

    if (WANTSP) {
      LWMIN = max(1, 2 * NN);
      LIWMIN = max(1, NN);
    } else if (lsame(JOB, 'N')) {
      LWMIN = max(1, N);
      LIWMIN = 1;
    } else if (lsame(JOB, 'E')) {
      LWMIN = max(1, NN);
      LIWMIN = 1;
    }

    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -15;
    } else if (LIWORK < LIWMIN && !LQUERY) {
      INFO.value = -17;
    }
  }

  if (INFO.value == 0) {
    WORK[1] = LWMIN.toDouble();
    IWORK[1] = LIWMIN;
  }

  if (INFO.value != 0) {
    xerbla('DTRSEN', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible.

  if (M.value == N || M.value == 0) {
    if (WANTS) S.value = ONE;
    if (WANTSP) SEP.value = dlange('1', N, N, T, LDT, WORK);
  } else {
    // Collect the selected blocks at the top-left corner of T.

    KS.value = 0;
    PAIR = false;
    var error = false;
    for (K = 1; K <= N; K++) {
      if (PAIR) {
        PAIR = false;
      } else {
        SWAP = SELECT[K];
        if (K < N) {
          if (T[K + 1][K] != ZERO) {
            PAIR = true;
            SWAP = SWAP || SELECT[K + 1];
          }
        }
        if (SWAP) {
          KS.value++;

          // Swap the K-th block to position KS.

          IERR.value = 0;
          KK.value = K;
          if (K != KS.value) {
            dtrexc(COMPQ, N, T, LDT, Q, LDQ, KK, KS, WORK, IERR);
          }
          if (IERR.value == 1 || IERR.value == 2) {
            // Blocks too close to swap: exit.

            INFO.value = 1;
            if (WANTS) S.value = ZERO;
            if (WANTSP) SEP.value = ZERO;
            error = true;
            break;
          }
          if (PAIR) KS.value++;
        }
      }
    }

    if (!error) {
      if (WANTS) {
        // Solve Sylvester equation for R:

        // T11*R - R*T22 = scale*T12

        dlacpy('F', N1, N2, T(1, N1 + 1), LDT, WORK.asMatrix(N1), N1);
        dtrsyl('N', 'N', -1, N1, N2, T, LDT, T(N1 + 1, N1 + 1), LDT,
            WORK.asMatrix(N1), N1, SCALE, IERR);

        // Estimate the reciprocal of the condition number of the cluster
        // of eigenvalues.

        RNORM = dlange('F', N1, N2, WORK.asMatrix(N1), N1, WORK);
        if (RNORM == ZERO) {
          S.value = ONE;
        } else {
          S.value = SCALE.value /
              (sqrt(SCALE.value * SCALE.value / RNORM + RNORM) * sqrt(RNORM));
        }
      }

      if (WANTSP) {
        // Estimate sep(T11,T22).

        EST.value = ZERO;
        KASE.value = 0;
        while (true) {
          dlacn2(NN, WORK(NN + 1), WORK, IWORK, EST, KASE, ISAVE);
          if (KASE.value == 0) break;
          if (KASE.value == 1) {
            // Solve  T11*R - R*T22 = scale*X.

            dtrsyl('N', 'N', -1, N1, N2, T, LDT, T(N1 + 1, N1 + 1), LDT,
                WORK.asMatrix(N1), N1, SCALE, IERR);
          } else {
            // Solve T11**T*R - R*T22**T = scale*X.

            dtrsyl('T', 'T', -1, N1, N2, T, LDT, T(N1 + 1, N1 + 1), LDT,
                WORK.asMatrix(N1), N1, SCALE, IERR);
          }
        }

        SEP.value = SCALE.value / EST.value;
      }
    }
  }

  // Store the output eigenvalues in WR and WI.
  for (K = 1; K <= N; K++) {
    WR[K] = T[K][K];
    WI[K] = ZERO;
  }
  for (K = 1; K <= N - 1; K++) {
    if (T[K + 1][K] != ZERO) {
      WI[K] = sqrt(T[K][K + 1].abs()) * sqrt(T[K + 1][K].abs());
      WI[K + 1] = -WI[K];
    }
  }

  WORK[1] = LWMIN.toDouble();
  IWORK[1] = LIWMIN;
}
