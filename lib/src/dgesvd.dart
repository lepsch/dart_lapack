// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dbdsqr.dart';
import 'package:dart_lapack/src/dgebrd.dart';
import 'package:dart_lapack/src/dgelqf.dart';
import 'package:dart_lapack/src/dgeqrf.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dlascl.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dorgbr.dart';
import 'package:dart_lapack/src/dorglq.dart';
import 'package:dart_lapack/src/dorgqr.dart';
import 'package:dart_lapack/src/dormbr.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dgesvd(
  final String JOBU,
  final String JOBVT,
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> S_,
  final Matrix<double> U_,
  final int LDU,
  final Matrix<double> VT_,
  final int LDVT,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final S = S_.having();
  final U = U_.having(ld: LDU);
  final VT = VT_.having(ld: LDVT);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool LQUERY,
      WNTUA,
      WNTUAS,
      WNTUN,
      WNTUO,
      WNTUS,
      WNTVA,
      WNTVAS,
      WNTVN,
      WNTVO,
      WNTVS;
  int BDSPAC = 0,
      BLK,
      CHUNK,
      I,
      IE = 0,
      IR,
      ISCL,
      ITAU,
      ITAUP,
      ITAUQ,
      IU,
      IWORK,
      LDWRKR,
      LDWRKU,
      MAXWRK = 0,
      MINMN,
      MINWRK,
      MNTHR = 0,
      NCU = 0,
      NCVT = 0,
      NRU = 0,
      NRVT = 0,
      WRKBL = 0;
  int LWORK_DGEQRF,
      LWORK_DORGQR_N,
      LWORK_DORGQR_M,
      LWORK_DGEBRD,
      LWORK_DORGBR_P,
      LWORK_DORGBR_Q,
      LWORK_DGELQF,
      LWORK_DORGLQ_N,
      LWORK_DORGLQ_M;
  double ANRM, BIGNUM, EPS, SMLNUM;
  final DUM = Array<double>(1);
  final IERR = Box(0);

  // Test the input arguments

  INFO.value = 0;
  MINMN = min(M, N);
  WNTUA = lsame(JOBU, 'A');
  WNTUS = lsame(JOBU, 'S');
  WNTUAS = WNTUA || WNTUS;
  WNTUO = lsame(JOBU, 'O');
  WNTUN = lsame(JOBU, 'N');
  WNTVA = lsame(JOBVT, 'A');
  WNTVS = lsame(JOBVT, 'S');
  WNTVAS = WNTVA || WNTVS;
  WNTVO = lsame(JOBVT, 'O');
  WNTVN = lsame(JOBVT, 'N');
  LQUERY = (LWORK == -1);

  if (!(WNTUA || WNTUS || WNTUO || WNTUN)) {
    INFO.value = -1;
  } else if (!(WNTVA || WNTVS || WNTVO || WNTVN) || (WNTVO && WNTUO)) {
    INFO.value = -2;
  } else if (M < 0) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (LDA < max(1, M)) {
    INFO.value = -6;
  } else if (LDU < 1 || (WNTUAS && LDU < M)) {
    INFO.value = -9;
  } else if (LDVT < 1 || (WNTVA && LDVT < N) || (WNTVS && LDVT < MINMN)) {
    INFO.value = -11;
  }

  // Compute workspace
  //  (Note: Comments in the code beginning "Workspace:" describe the
  //   minimal amount of workspace needed at that point in the code,
  //   as well as the preferred amount for good performance.
  //   NB refers to the optimal block size for the immediately
  //   following subroutine, as returned by ILAENV.)

  if (INFO.value == 0) {
    MINWRK = 1;
    MAXWRK = 1;
    if (M >= N && MINMN > 0) {
      // Compute space needed for DBDSQR

      MNTHR = ilaenv(6, 'DGESVD', JOBU + JOBVT, M, N, 0, 0);
      BDSPAC = 5 * N;
      // Compute space needed for DGEQRF
      dgeqrf(M, N, A, LDA, DUM(1), DUM(1), -1, IERR);
      LWORK_DGEQRF = DUM[1].toInt();
      // Compute space needed for DORGQR
      dorgqr(M, N, N, A, LDA, DUM(1), DUM(1), -1, IERR);
      LWORK_DORGQR_N = DUM[1].toInt();
      dorgqr(M, M, N, A, LDA, DUM(1), DUM(1), -1, IERR);
      LWORK_DORGQR_M = DUM[1].toInt();
      // Compute space needed for DGEBRD
      dgebrd(N, N, A, LDA, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR);
      LWORK_DGEBRD = DUM[1].toInt();
      // Compute space needed for DORGBR P
      dorgbr('P', N, N, N, A, LDA, DUM(1), DUM(1), -1, IERR);
      LWORK_DORGBR_P = DUM[1].toInt();
      // Compute space needed for DORGBR Q
      dorgbr('Q', N, N, N, A, LDA, DUM(1), DUM(1), -1, IERR);
      LWORK_DORGBR_Q = DUM[1].toInt();

      if (M >= MNTHR) {
        if (WNTUN) {
          // Path 1 (M much larger than N, JOBU='N')

          MAXWRK = N + LWORK_DGEQRF;
          MAXWRK = max(MAXWRK, 3 * N + LWORK_DGEBRD);
          if (WNTVO || WNTVAS) MAXWRK = max(MAXWRK, 3 * N + LWORK_DORGBR_P);
          MAXWRK = max(MAXWRK, BDSPAC);
          MINWRK = max(4 * N, BDSPAC);
        } else if (WNTUO && WNTVN) {
          // Path 2 (M much larger than N, JOBU='O', JOBVT='N')

          WRKBL = N + LWORK_DGEQRF;
          WRKBL = max(WRKBL, N + LWORK_DORGQR_N);
          WRKBL = max(WRKBL, 3 * N + LWORK_DGEBRD);
          WRKBL = max(WRKBL, 3 * N + LWORK_DORGBR_Q);
          WRKBL = max(WRKBL, BDSPAC);
          MAXWRK = max(N * N + WRKBL, N * N + M * N + N);
          MINWRK = max(3 * N + M, BDSPAC);
        } else if (WNTUO && WNTVAS) {
          // Path 3 (M much larger than N, JOBU='O', JOBVT='S' or
          // 'A')

          WRKBL = N + LWORK_DGEQRF;
          WRKBL = max(WRKBL, N + LWORK_DORGQR_N);
          WRKBL = max(WRKBL, 3 * N + LWORK_DGEBRD);
          WRKBL = max(WRKBL, 3 * N + LWORK_DORGBR_Q);
          WRKBL = max(WRKBL, 3 * N + LWORK_DORGBR_P);
          WRKBL = max(WRKBL, BDSPAC);
          MAXWRK = max(N * N + WRKBL, N * N + M * N + N);
          MINWRK = max(3 * N + M, BDSPAC);
        } else if (WNTUS && WNTVN) {
          // Path 4 (M much larger than N, JOBU='S', JOBVT='N')

          WRKBL = N + LWORK_DGEQRF;
          WRKBL = max(WRKBL, N + LWORK_DORGQR_N);
          WRKBL = max(WRKBL, 3 * N + LWORK_DGEBRD);
          WRKBL = max(WRKBL, 3 * N + LWORK_DORGBR_Q);
          WRKBL = max(WRKBL, BDSPAC);
          MAXWRK = N * N + WRKBL;
          MINWRK = max(3 * N + M, BDSPAC);
        } else if (WNTUS && WNTVO) {
          // Path 5 (M much larger than N, JOBU='S', JOBVT='O')

          WRKBL = N + LWORK_DGEQRF;
          WRKBL = max(WRKBL, N + LWORK_DORGQR_N);
          WRKBL = max(WRKBL, 3 * N + LWORK_DGEBRD);
          WRKBL = max(WRKBL, 3 * N + LWORK_DORGBR_Q);
          WRKBL = max(WRKBL, 3 * N + LWORK_DORGBR_P);
          WRKBL = max(WRKBL, BDSPAC);
          MAXWRK = 2 * N * N + WRKBL;
          MINWRK = max(3 * N + M, BDSPAC);
        } else if (WNTUS && WNTVAS) {
          // Path 6 (M much larger than N, JOBU='S', JOBVT='S' or
          // 'A')

          WRKBL = N + LWORK_DGEQRF;
          WRKBL = max(WRKBL, N + LWORK_DORGQR_N);
          WRKBL = max(WRKBL, 3 * N + LWORK_DGEBRD);
          WRKBL = max(WRKBL, 3 * N + LWORK_DORGBR_Q);
          WRKBL = max(WRKBL, 3 * N + LWORK_DORGBR_P);
          WRKBL = max(WRKBL, BDSPAC);
          MAXWRK = N * N + WRKBL;
          MINWRK = max(3 * N + M, BDSPAC);
        } else if (WNTUA && WNTVN) {
          // Path 7 (M much larger than N, JOBU='A', JOBVT='N')

          WRKBL = N + LWORK_DGEQRF;
          WRKBL = max(WRKBL, N + LWORK_DORGQR_M);
          WRKBL = max(WRKBL, 3 * N + LWORK_DGEBRD);
          WRKBL = max(WRKBL, 3 * N + LWORK_DORGBR_Q);
          WRKBL = max(WRKBL, BDSPAC);
          MAXWRK = N * N + WRKBL;
          MINWRK = max(3 * N + M, BDSPAC);
        } else if (WNTUA && WNTVO) {
          // Path 8 (M much larger than N, JOBU='A', JOBVT='O')

          WRKBL = N + LWORK_DGEQRF;
          WRKBL = max(WRKBL, N + LWORK_DORGQR_M);
          WRKBL = max(WRKBL, 3 * N + LWORK_DGEBRD);
          WRKBL = max(WRKBL, 3 * N + LWORK_DORGBR_Q);
          WRKBL = max(WRKBL, 3 * N + LWORK_DORGBR_P);
          WRKBL = max(WRKBL, BDSPAC);
          MAXWRK = 2 * N * N + WRKBL;
          MINWRK = max(3 * N + M, BDSPAC);
        } else if (WNTUA && WNTVAS) {
          // Path 9 (M much larger than N, JOBU='A', JOBVT='S' or
          // 'A')

          WRKBL = N + LWORK_DGEQRF;
          WRKBL = max(WRKBL, N + LWORK_DORGQR_M);
          WRKBL = max(WRKBL, 3 * N + LWORK_DGEBRD);
          WRKBL = max(WRKBL, 3 * N + LWORK_DORGBR_Q);
          WRKBL = max(WRKBL, 3 * N + LWORK_DORGBR_P);
          WRKBL = max(WRKBL, BDSPAC);
          MAXWRK = N * N + WRKBL;
          MINWRK = max(3 * N + M, BDSPAC);
        }
      } else {
        // Path 10 (M at least N, but not much larger)

        dgebrd(M, N, A, LDA, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR);
        LWORK_DGEBRD = DUM[1].toInt();
        MAXWRK = 3 * N + LWORK_DGEBRD;
        if (WNTUS || WNTUO) {
          dorgbr('Q', M, N, N, A, LDA, DUM(1), DUM(1), -1, IERR);
          LWORK_DORGBR_Q = DUM[1].toInt();
          MAXWRK = max(MAXWRK, 3 * N + LWORK_DORGBR_Q);
        }
        if (WNTUA) {
          dorgbr('Q', M, M, N, A, LDA, DUM(1), DUM(1), -1, IERR);
          LWORK_DORGBR_Q = DUM[1].toInt();
          MAXWRK = max(MAXWRK, 3 * N + LWORK_DORGBR_Q);
        }
        if (!WNTVN) {
          MAXWRK = max(MAXWRK, 3 * N + LWORK_DORGBR_P);
        }
        MAXWRK = max(MAXWRK, BDSPAC);
        MINWRK = max(3 * N + M, BDSPAC);
      }
    } else if (MINMN > 0) {
      // Compute space needed for DBDSQR

      MNTHR = ilaenv(6, 'DGESVD', JOBU + JOBVT, M, N, 0, 0);
      BDSPAC = 5 * M;
      // Compute space needed for DGELQF
      dgelqf(M, N, A, LDA, DUM(1), DUM(1), -1, IERR);
      LWORK_DGELQF = DUM[1].toInt();
      // Compute space needed for DORGLQ
      dorglq(N, N, M, DUM(1).asMatrix(N), N, DUM(1), DUM(1), -1, IERR);
      LWORK_DORGLQ_N = DUM[1].toInt();
      dorglq(M, N, M, A, LDA, DUM(1), DUM(1), -1, IERR);
      LWORK_DORGLQ_M = DUM[1].toInt();
      // Compute space needed for DGEBRD
      dgebrd(M, M, A, LDA, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR);
      LWORK_DGEBRD = DUM[1].toInt();
      // Compute space needed for DORGBR P
      dorgbr('P', M, M, M, A, N, DUM(1), DUM(1), -1, IERR);
      LWORK_DORGBR_P = DUM[1].toInt();
      // Compute space needed for DORGBR Q
      dorgbr('Q', M, M, M, A, N, DUM(1), DUM(1), -1, IERR);
      LWORK_DORGBR_Q = DUM[1].toInt();
      if (N >= MNTHR) {
        if (WNTVN) {
          // Path 1t(N much larger than M, JOBVT='N')

          MAXWRK = M + LWORK_DGELQF;
          MAXWRK = max(MAXWRK, 3 * M + LWORK_DGEBRD);
          if (WNTUO || WNTUAS) MAXWRK = max(MAXWRK, 3 * M + LWORK_DORGBR_Q);
          MAXWRK = max(MAXWRK, BDSPAC);
          MINWRK = max(4 * M, BDSPAC);
        } else if (WNTVO && WNTUN) {
          // Path 2t(N much larger than M, JOBU='N', JOBVT='O')

          WRKBL = M + LWORK_DGELQF;
          WRKBL = max(WRKBL, M + LWORK_DORGLQ_M);
          WRKBL = max(WRKBL, 3 * M + LWORK_DGEBRD);
          WRKBL = max(WRKBL, 3 * M + LWORK_DORGBR_P);
          WRKBL = max(WRKBL, BDSPAC);
          MAXWRK = max(M * M + WRKBL, M * M + M * N + M);
          MINWRK = max(3 * M + N, BDSPAC);
        } else if (WNTVO && WNTUAS) {
          // Path 3t(N much larger than M, JOBU='S' or 'A',
          // JOBVT='O')

          WRKBL = M + LWORK_DGELQF;
          WRKBL = max(WRKBL, M + LWORK_DORGLQ_M);
          WRKBL = max(WRKBL, 3 * M + LWORK_DGEBRD);
          WRKBL = max(WRKBL, 3 * M + LWORK_DORGBR_P);
          WRKBL = max(WRKBL, 3 * M + LWORK_DORGBR_Q);
          WRKBL = max(WRKBL, BDSPAC);
          MAXWRK = max(M * M + WRKBL, M * M + M * N + M);
          MINWRK = max(3 * M + N, BDSPAC);
        } else if (WNTVS && WNTUN) {
          // Path 4t(N much larger than M, JOBU='N', JOBVT='S')

          WRKBL = M + LWORK_DGELQF;
          WRKBL = max(WRKBL, M + LWORK_DORGLQ_M);
          WRKBL = max(WRKBL, 3 * M + LWORK_DGEBRD);
          WRKBL = max(WRKBL, 3 * M + LWORK_DORGBR_P);
          WRKBL = max(WRKBL, BDSPAC);
          MAXWRK = M * M + WRKBL;
          MINWRK = max(3 * M + N, BDSPAC);
        } else if (WNTVS && WNTUO) {
          // Path 5t(N much larger than M, JOBU='O', JOBVT='S')

          WRKBL = M + LWORK_DGELQF;
          WRKBL = max(WRKBL, M + LWORK_DORGLQ_M);
          WRKBL = max(WRKBL, 3 * M + LWORK_DGEBRD);
          WRKBL = max(WRKBL, 3 * M + LWORK_DORGBR_P);
          WRKBL = max(WRKBL, 3 * M + LWORK_DORGBR_Q);
          WRKBL = max(WRKBL, BDSPAC);
          MAXWRK = 2 * M * M + WRKBL;
          MINWRK = max(3 * M + N, BDSPAC);
        } else if (WNTVS && WNTUAS) {
          // Path 6t(N much larger than M, JOBU='S' or 'A',
          // JOBVT='S')

          WRKBL = M + LWORK_DGELQF;
          WRKBL = max(WRKBL, M + LWORK_DORGLQ_M);
          WRKBL = max(WRKBL, 3 * M + LWORK_DGEBRD);
          WRKBL = max(WRKBL, 3 * M + LWORK_DORGBR_P);
          WRKBL = max(WRKBL, 3 * M + LWORK_DORGBR_Q);
          WRKBL = max(WRKBL, BDSPAC);
          MAXWRK = M * M + WRKBL;
          MINWRK = max(3 * M + N, BDSPAC);
        } else if (WNTVA && WNTUN) {
          // Path 7t(N much larger than M, JOBU='N', JOBVT='A')

          WRKBL = M + LWORK_DGELQF;
          WRKBL = max(WRKBL, M + LWORK_DORGLQ_N);
          WRKBL = max(WRKBL, 3 * M + LWORK_DGEBRD);
          WRKBL = max(WRKBL, 3 * M + LWORK_DORGBR_P);
          WRKBL = max(WRKBL, BDSPAC);
          MAXWRK = M * M + WRKBL;
          MINWRK = max(3 * M + N, BDSPAC);
        } else if (WNTVA && WNTUO) {
          // Path 8t(N much larger than M, JOBU='O', JOBVT='A')

          WRKBL = M + LWORK_DGELQF;
          WRKBL = max(WRKBL, M + LWORK_DORGLQ_N);
          WRKBL = max(WRKBL, 3 * M + LWORK_DGEBRD);
          WRKBL = max(WRKBL, 3 * M + LWORK_DORGBR_P);
          WRKBL = max(WRKBL, 3 * M + LWORK_DORGBR_Q);
          WRKBL = max(WRKBL, BDSPAC);
          MAXWRK = 2 * M * M + WRKBL;
          MINWRK = max(3 * M + N, BDSPAC);
        } else if (WNTVA && WNTUAS) {
          // Path 9t(N much larger than M, JOBU='S' or 'A',
          // JOBVT='A')

          WRKBL = M + LWORK_DGELQF;
          WRKBL = max(WRKBL, M + LWORK_DORGLQ_N);
          WRKBL = max(WRKBL, 3 * M + LWORK_DGEBRD);
          WRKBL = max(WRKBL, 3 * M + LWORK_DORGBR_P);
          WRKBL = max(WRKBL, 3 * M + LWORK_DORGBR_Q);
          WRKBL = max(WRKBL, BDSPAC);
          MAXWRK = M * M + WRKBL;
          MINWRK = max(3 * M + N, BDSPAC);
        }
      } else {
        // Path 10t(N greater than M, but not much larger)

        dgebrd(M, N, A, LDA, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR);
        LWORK_DGEBRD = DUM[1].toInt();
        MAXWRK = 3 * M + LWORK_DGEBRD;
        if (WNTVS || WNTVO) {
          // Compute space needed for DORGBR P
          dorgbr('P', M, N, M, A, N, DUM(1), DUM(1), -1, IERR);
          LWORK_DORGBR_P = DUM[1].toInt();
          MAXWRK = max(MAXWRK, 3 * M + LWORK_DORGBR_P);
        }
        if (WNTVA) {
          dorgbr('P', N, N, M, A, N, DUM(1), DUM(1), -1, IERR);
          LWORK_DORGBR_P = DUM[1].toInt();
          MAXWRK = max(MAXWRK, 3 * M + LWORK_DORGBR_P);
        }
        if (!WNTUN) {
          MAXWRK = max(MAXWRK, 3 * M + LWORK_DORGBR_Q);
        }
        MAXWRK = max(MAXWRK, BDSPAC);
        MINWRK = max(3 * M + N, BDSPAC);
      }
    }
    MAXWRK = max(MAXWRK, MINWRK);
    WORK[1] = MAXWRK.toDouble();

    if (LWORK < MINWRK && !LQUERY) {
      INFO.value = -13;
    }
  }

  if (INFO.value != 0) {
    xerbla('DGESVD', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) {
    return;
  }

  // Get machine constants

  EPS = dlamch('P');
  SMLNUM = sqrt(dlamch('S')) / EPS;
  BIGNUM = ONE / SMLNUM;

  // Scale A if max element outside range [SMLNUM,BIGNUM]

  ANRM = dlange('M', M, N, A, LDA, DUM);
  ISCL = 0;
  if (ANRM > ZERO && ANRM < SMLNUM) {
    ISCL = 1;
    dlascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, IERR);
  } else if (ANRM > BIGNUM) {
    ISCL = 1;
    dlascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, IERR);
  }

  if (M >= N) {
    // A has at least as many rows as columns. If A has sufficiently
    // more rows than columns, first reduce using the QR
    // decomposition (if sufficient workspace available)

    if (M >= MNTHR) {
      if (WNTUN) {
        // Path 1 (M much larger than N, JOBU='N')
        // No left singular vectors to be computed

        ITAU = 1;
        IWORK = ITAU + N;

        // Compute A=Q*R
        // (Workspace: need 2*N, prefer N + N*NB)

        dgeqrf(M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

        // Zero out below R

        if (N > 1) {
          dlaset('L', N - 1, N - 1, ZERO, ZERO, A(2, 1), LDA);
        }
        IE = 1;
        ITAUQ = IE + N;
        ITAUP = ITAUQ + N;
        IWORK = ITAUP + N;

        // Bidiagonalize R in A
        // (Workspace: need 4*N, prefer 3*N + 2*N*NB)

        dgebrd(N, N, A, LDA, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP), WORK(IWORK),
            LWORK - IWORK + 1, IERR);
        NCVT = 0;
        if (WNTVO || WNTVAS) {
          // If right singular vectors desired, generate P'.
          // (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)

          dorgbr('P', N, N, N, A, LDA, WORK(ITAUP), WORK(IWORK),
              LWORK - IWORK + 1, IERR);
          NCVT = N;
        }
        IWORK = IE + N;

        // Perform bidiagonal QR iteration, computing right
        // singular vectors of A in A if desired
        // (Workspace: need BDSPAC)

        dbdsqr('U', N, NCVT, 0, 0, S, WORK(IE), A, LDA, DUM.asMatrix(1), 1,
            DUM.asMatrix(1), 1, WORK(IWORK), INFO);

        // If right singular vectors desired in VT, copy them there

        if (WNTVAS) dlacpy('F', N, N, A, LDA, VT, LDVT);
      } else if (WNTUO && WNTVN) {
        // Path 2 (M much larger than N, JOBU='O', JOBVT='N')
        // N left singular vectors to be overwritten on A and
        // no right singular vectors to be computed

        if (LWORK >= N * N + max(4 * N, BDSPAC)) {
          // Sufficient workspace for a fast algorithm

          IR = 1;
          if (LWORK >= max(WRKBL, LDA * N + N) + LDA * N) {
            // WORK(IU) is LDA by N, WORK(IR) is LDA by N

            LDWRKU = LDA;
            LDWRKR = LDA;
          } else if (LWORK >= max(WRKBL, LDA * N + N) + N * N) {
            // WORK(IU) is LDA by N, WORK(IR) is N by N

            LDWRKU = LDA;
            LDWRKR = N;
          } else {
            // WORK(IU) is LDWRKU by N, WORK(IR) is N by N

            LDWRKU = (LWORK - N * N - N) ~/ N;
            LDWRKR = N;
          }
          ITAU = IR + LDWRKR * N;
          IWORK = ITAU + N;

          // Compute A=Q*R
          // (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)

          dgeqrf(
              M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Copy R to WORK(IR) and zero out below it

          dlacpy('U', N, N, A, LDA, WORK(IR).asMatrix(LDWRKR), LDWRKR);
          dlaset('L', N - 1, N - 1, ZERO, ZERO, WORK(IR + 1).asMatrix(LDWRKR),
              LDWRKR);

          // Generate Q in A
          // (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)

          dorgqr(M, N, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
              IERR);
          IE = ITAU;
          ITAUQ = IE + N;
          ITAUP = ITAUQ + N;
          IWORK = ITAUP + N;

          // Bidiagonalize R in WORK(IR)
          // (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB)

          dgebrd(N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, S, WORK(IE),
              WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Generate left vectors bidiagonalizing R
          // (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB)

          dorgbr('Q', N, N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUQ),
              WORK(IWORK), LWORK - IWORK + 1, IERR);
          IWORK = IE + N;

          // Perform bidiagonal QR iteration, computing left
          // singular vectors of R in WORK(IR)
          // (Workspace: need N*N + BDSPAC)

          dbdsqr(
              'U',
              N,
              0,
              N,
              0,
              S,
              WORK(IE),
              DUM.asMatrix(1),
              1,
              WORK(IR).asMatrix(LDWRKR),
              LDWRKR,
              DUM.asMatrix(1),
              1,
              WORK(IWORK),
              INFO);
          IU = IE + N;

          // Multiply Q in A by left singular vectors of R in
          // WORK(IR), storing result in WORK(IU) and copying to A
          // (Workspace: need N*N + 2*N, prefer N*N + M*N + N)

          for (I = 1; I <= M; I += LDWRKU) {
            CHUNK = min(M - I + 1, LDWRKU);
            dgemm(
                'N',
                'N',
                CHUNK,
                N,
                N,
                ONE,
                A(I, 1),
                LDA,
                WORK(IR).asMatrix(LDWRKR),
                LDWRKR,
                ZERO,
                WORK(IU).asMatrix(LDWRKU),
                LDWRKU);
            dlacpy(
                'F', CHUNK, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, A(I, 1), LDA);
          }
        } else {
          // Insufficient workspace for a fast algorithm

          IE = 1;
          ITAUQ = IE + N;
          ITAUP = ITAUQ + N;
          IWORK = ITAUP + N;

          // Bidiagonalize A
          // (Workspace: need 3*N + M, prefer 3*N + (M + N)*NB)

          dgebrd(M, N, A, LDA, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP),
              WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Generate left vectors bidiagonalizing A
          // (Workspace: need 4*N, prefer 3*N + N*NB)

          dorgbr('Q', M, N, N, A, LDA, WORK(ITAUQ), WORK(IWORK),
              LWORK - IWORK + 1, IERR);
          IWORK = IE + N;

          // Perform bidiagonal QR iteration, computing left
          // singular vectors of A in A
          // (Workspace: need BDSPAC)

          dbdsqr('U', N, 0, M, 0, S, WORK(IE), DUM.asMatrix(1), 1, A, LDA,
              DUM.asMatrix(1), 1, WORK(IWORK), INFO);
        }
      } else if (WNTUO && WNTVAS) {
        // Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A')
        // N left singular vectors to be overwritten on A and
        // N right singular vectors to be computed in VT

        if (LWORK >= N * N + max(4 * N, BDSPAC)) {
          // Sufficient workspace for a fast algorithm

          IR = 1;
          if (LWORK >= max(WRKBL, LDA * N + N) + LDA * N) {
            // WORK(IU) is LDA by N and WORK(IR) is LDA by N

            LDWRKU = LDA;
            LDWRKR = LDA;
          } else if (LWORK >= max(WRKBL, LDA * N + N) + N * N) {
            // WORK(IU) is LDA by N and WORK(IR) is N by N

            LDWRKU = LDA;
            LDWRKR = N;
          } else {
            // WORK(IU) is LDWRKU by N and WORK(IR) is N by N

            LDWRKU = (LWORK - N * N - N) ~/ N;
            LDWRKR = N;
          }
          ITAU = IR + LDWRKR * N;
          IWORK = ITAU + N;

          // Compute A=Q*R
          // (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)

          dgeqrf(
              M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Copy R to VT, zeroing out below it

          dlacpy('U', N, N, A, LDA, VT, LDVT);
          if (N > 1) dlaset('L', N - 1, N - 1, ZERO, ZERO, VT(2, 1), LDVT);

          // Generate Q in A
          // (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)

          dorgqr(M, N, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
              IERR);
          IE = ITAU;
          ITAUQ = IE + N;
          ITAUP = ITAUQ + N;
          IWORK = ITAUP + N;

          // Bidiagonalize R in VT, copying result to WORK(IR)
          // (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB)

          dgebrd(N, N, VT, LDVT, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP),
              WORK(IWORK), LWORK - IWORK + 1, IERR);
          dlacpy('L', N, N, VT, LDVT, WORK(IR).asMatrix(LDWRKR), LDWRKR);

          // Generate left vectors bidiagonalizing R in WORK(IR)
          // (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB)

          dorgbr('Q', N, N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUQ),
              WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Generate right vectors bidiagonalizing R in VT
          // (Workspace: need N*N + 4*N-1, prefer N*N + 3*N + (N-1)*NB)

          dorgbr('P', N, N, N, VT, LDVT, WORK(ITAUP), WORK(IWORK),
              LWORK - IWORK + 1, IERR);
          IWORK = IE + N;

          // Perform bidiagonal QR iteration, computing left
          // singular vectors of R in WORK(IR) and computing right
          // singular vectors of R in VT
          // (Workspace: need N*N + BDSPAC)

          dbdsqr(
              'U',
              N,
              N,
              N,
              0,
              S,
              WORK(IE),
              VT,
              LDVT,
              WORK(IR).asMatrix(LDWRKR),
              LDWRKR,
              DUM.asMatrix(1),
              1,
              WORK(IWORK),
              INFO);
          IU = IE + N;

          // Multiply Q in A by left singular vectors of R in
          // WORK(IR), storing result in WORK(IU) and copying to A
          // (Workspace: need N*N + 2*N, prefer N*N + M*N + N)

          for (I = 1; I <= M; I += LDWRKU) {
            CHUNK = min(M - I + 1, LDWRKU);
            dgemm(
                'N',
                'N',
                CHUNK,
                N,
                N,
                ONE,
                A(I, 1),
                LDA,
                WORK(IR).asMatrix(LDWRKR),
                LDWRKR,
                ZERO,
                WORK(IU).asMatrix(LDWRKU),
                LDWRKU);
            dlacpy(
                'F', CHUNK, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, A(I, 1), LDA);
          }
        } else {
          // Insufficient workspace for a fast algorithm

          ITAU = 1;
          IWORK = ITAU + N;

          // Compute A=Q*R
          // (Workspace: need 2*N, prefer N + N*NB)

          dgeqrf(
              M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Copy R to VT, zeroing out below it

          dlacpy('U', N, N, A, LDA, VT, LDVT);
          if (N > 1) dlaset('L', N - 1, N - 1, ZERO, ZERO, VT(2, 1), LDVT);

          // Generate Q in A
          // (Workspace: need 2*N, prefer N + N*NB)

          dorgqr(M, N, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
              IERR);
          IE = ITAU;
          ITAUQ = IE + N;
          ITAUP = ITAUQ + N;
          IWORK = ITAUP + N;

          // Bidiagonalize R in VT
          // (Workspace: need 4*N, prefer 3*N + 2*N*NB)

          dgebrd(N, N, VT, LDVT, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP),
              WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Multiply Q in A by left vectors bidiagonalizing R
          // (Workspace: need 3*N + M, prefer 3*N + M*NB)

          dormbr('Q', 'R', 'N', M, N, N, VT, LDVT, WORK(ITAUQ), A, LDA,
              WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Generate right vectors bidiagonalizing R in VT
          // (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)

          dorgbr('P', N, N, N, VT, LDVT, WORK(ITAUP), WORK(IWORK),
              LWORK - IWORK + 1, IERR);
          IWORK = IE + N;

          // Perform bidiagonal QR iteration, computing left
          // singular vectors of A in A and computing right
          // singular vectors of A in VT
          // (Workspace: need BDSPAC)

          dbdsqr('U', N, N, M, 0, S, WORK(IE), VT, LDVT, A, LDA,
              DUM.asMatrix(1), 1, WORK(IWORK), INFO);
        }
      } else if (WNTUS) {
        if (WNTVN) {
          // Path 4 (M much larger than N, JOBU='S', JOBVT='N')
          // N left singular vectors to be computed in U and
          // no right singular vectors to be computed

          if (LWORK >= N * N + max(4 * N, BDSPAC)) {
            // Sufficient workspace for a fast algorithm

            IR = 1;
            if (LWORK >= WRKBL + LDA * N) {
              // WORK(IR) is LDA by N

              LDWRKR = LDA;
            } else {
              // WORK(IR) is N by N

              LDWRKR = N;
            }
            ITAU = IR + LDWRKR * N;
            IWORK = ITAU + N;

            // Compute A=Q*R
            // (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)

            dgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Copy R to WORK(IR), zeroing out below it

            dlacpy('U', N, N, A, LDA, WORK(IR).asMatrix(LDWRKR), LDWRKR);
            dlaset('L', N - 1, N - 1, ZERO, ZERO, WORK(IR + 1).asMatrix(LDWRKR),
                LDWRKR);

            // Generate Q in A
            // (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)

            dorgqr(M, N, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);
            IE = ITAU;
            ITAUQ = IE + N;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Bidiagonalize R in WORK(IR)
            // (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB)

            dgebrd(N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, S, WORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate left vectors bidiagonalizing R in WORK(IR)
            // (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB)

            dorgbr('Q', N, N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUQ),
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of R in WORK(IR)
            // (Workspace: need N*N + BDSPAC)

            dbdsqr(
                'U',
                N,
                0,
                N,
                0,
                S,
                WORK(IE),
                DUM.asMatrix(1),
                1,
                WORK(IR).asMatrix(LDWRKR),
                LDWRKR,
                DUM.asMatrix(1),
                1,
                WORK(IWORK),
                INFO);

            // Multiply Q in A by left singular vectors of R in
            // WORK(IR), storing result in U
            // (Workspace: need N*N)

            dgemm('N', 'N', M, N, N, ONE, A, LDA, WORK(IR).asMatrix(LDWRKR),
                LDWRKR, ZERO, U, LDU);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + N;

            // Compute A=Q*R, copying result to U
            // (Workspace: need 2*N, prefer N + N*NB)

            dgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('L', M, N, A, LDA, U, LDU);

            // Generate Q in U
            // (Workspace: need 2*N, prefer N + N*NB)

            dorgqr(M, N, N, U, LDU, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);
            IE = ITAU;
            ITAUQ = IE + N;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Zero out below R in A

            if (N > 1) {
              dlaset('L', N - 1, N - 1, ZERO, ZERO, A(2, 1), LDA);
            }

            // Bidiagonalize R in A
            // (Workspace: need 4*N, prefer 3*N + 2*N*NB)

            dgebrd(N, N, A, LDA, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply Q in U by left vectors bidiagonalizing R
            // (Workspace: need 3*N + M, prefer 3*N + M*NB)

            dormbr('Q', 'R', 'N', M, N, N, A, LDA, WORK(ITAUQ), U, LDU,
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of A in U
            // (Workspace: need BDSPAC)

            dbdsqr('U', N, 0, M, 0, S, WORK(IE), DUM.asMatrix(1), 1, U, LDU,
                DUM.asMatrix(1), 1, WORK(IWORK), INFO);
          }
        } else if (WNTVO) {
          // Path 5 (M much larger than N, JOBU='S', JOBVT='O')
          // N left singular vectors to be computed in U and
          // N right singular vectors to be overwritten on A

          if (LWORK >= 2 * N * N + max(4 * N, BDSPAC)) {
            // Sufficient workspace for a fast algorithm

            IU = 1;
            if (LWORK >= WRKBL + 2 * LDA * N) {
              // WORK(IU) is LDA by N and WORK(IR) is LDA by N

              LDWRKU = LDA;
              IR = IU + LDWRKU * N;
              LDWRKR = LDA;
            } else if (LWORK >= WRKBL + (LDA + N) * N) {
              // WORK(IU) is LDA by N and WORK(IR) is N by N

              LDWRKU = LDA;
              IR = IU + LDWRKU * N;
              LDWRKR = N;
            } else {
              // WORK(IU) is N by N and WORK(IR) is N by N

              LDWRKU = N;
              IR = IU + LDWRKU * N;
              LDWRKR = N;
            }
            ITAU = IR + LDWRKR * N;
            IWORK = ITAU + N;

            // Compute A=Q*R
            // (Workspace: need 2*N*N + 2*N, prefer 2*N*N + N + N*NB)

            dgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Copy R to WORK(IU), zeroing out below it

            dlacpy('U', N, N, A, LDA, WORK(IU).asMatrix(LDWRKU), LDWRKU);
            dlaset('L', N - 1, N - 1, ZERO, ZERO, WORK(IU + 1).asMatrix(LDWRKU),
                LDWRKU);

            // Generate Q in A
            // (Workspace: need 2*N*N + 2*N, prefer 2*N*N + N + N*NB)

            dorgqr(M, N, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);
            IE = ITAU;
            ITAUQ = IE + N;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Bidiagonalize R in WORK(IU), copying result to
            // WORK(IR)
            // (Workspace: need 2*N*N + 4*N,
            //             prefer 2*N*N+3*N+2*N*NB)

            dgebrd(N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, S, WORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('U', N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU,
                WORK(IR).asMatrix(LDWRKR), LDWRKR);

            // Generate left bidiagonalizing vectors in WORK(IU)
            // (Workspace: need 2*N*N + 4*N, prefer 2*N*N + 3*N + N*NB)

            dorgbr('Q', N, N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, WORK(ITAUQ),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate right bidiagonalizing vectors in WORK(IR)
            // (Workspace: need 2*N*N + 4*N-1,
            //             prefer 2*N*N+3*N+(N-1)*NB)

            dorgbr('P', N, N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of R in WORK(IU) and computing
            // right singular vectors of R in WORK(IR)
            // (Workspace: need 2*N*N + BDSPAC)

            dbdsqr(
                'U',
                N,
                N,
                N,
                0,
                S,
                WORK(IE),
                WORK(IR).asMatrix(LDWRKR),
                LDWRKR,
                WORK(IU).asMatrix(LDWRKU),
                LDWRKU,
                DUM.asMatrix(1),
                1,
                WORK(IWORK),
                INFO);

            // Multiply Q in A by left singular vectors of R in
            // WORK(IU), storing result in U
            // (Workspace: need N*N)

            dgemm('N', 'N', M, N, N, ONE, A, LDA, WORK(IU).asMatrix(LDWRKU),
                LDWRKU, ZERO, U, LDU);

            // Copy right singular vectors of R to A
            // (Workspace: need N*N)

            dlacpy('F', N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, A, LDA);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + N;

            // Compute A=Q*R, copying result to U
            // (Workspace: need 2*N, prefer N + N*NB)

            dgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('L', M, N, A, LDA, U, LDU);

            // Generate Q in U
            // (Workspace: need 2*N, prefer N + N*NB)

            dorgqr(M, N, N, U, LDU, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);
            IE = ITAU;
            ITAUQ = IE + N;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Zero out below R in A

            if (N > 1) {
              dlaset('L', N - 1, N - 1, ZERO, ZERO, A(2, 1), LDA);
            }

            // Bidiagonalize R in A
            // (Workspace: need 4*N, prefer 3*N + 2*N*NB)

            dgebrd(N, N, A, LDA, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply Q in U by left vectors bidiagonalizing R
            // (Workspace: need 3*N + M, prefer 3*N + M*NB)

            dormbr('Q', 'R', 'N', M, N, N, A, LDA, WORK(ITAUQ), U, LDU,
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate right vectors bidiagonalizing R in A
            // (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)

            dorgbr('P', N, N, N, A, LDA, WORK(ITAUP), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of A in U and computing right
            // singular vectors of A in A
            // (Workspace: need BDSPAC)

            dbdsqr('U', N, N, M, 0, S, WORK(IE), A, LDA, U, LDU,
                DUM.asMatrix(1), 1, WORK(IWORK), INFO);
          }
        } else if (WNTVAS) {
          // Path 6 (M much larger than N, JOBU='S', JOBVT='S'
          //         or 'A')
          // N left singular vectors to be computed in U and
          // N right singular vectors to be computed in VT

          if (LWORK >= N * N + max(4 * N, BDSPAC)) {
            // Sufficient workspace for a fast algorithm

            IU = 1;
            if (LWORK >= WRKBL + LDA * N) {
              // WORK(IU) is LDA by N

              LDWRKU = LDA;
            } else {
              // WORK(IU) is N by N

              LDWRKU = N;
            }
            ITAU = IU + LDWRKU * N;
            IWORK = ITAU + N;

            // Compute A=Q*R
            // (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)

            dgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Copy R to WORK(IU), zeroing out below it

            dlacpy('U', N, N, A, LDA, WORK(IU).asMatrix(LDWRKU), LDWRKU);
            dlaset('L', N - 1, N - 1, ZERO, ZERO, WORK(IU + 1).asMatrix(LDWRKU),
                LDWRKU);

            // Generate Q in A
            // (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)

            dorgqr(M, N, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);
            IE = ITAU;
            ITAUQ = IE + N;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Bidiagonalize R in WORK(IU), copying result to VT
            // (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB)

            dgebrd(N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, S, WORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('U', N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, VT, LDVT);

            // Generate left bidiagonalizing vectors in WORK(IU)
            // (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB)

            dorgbr('Q', N, N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, WORK(ITAUQ),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate right bidiagonalizing vectors in VT
            // (Workspace: need N*N + 4*N-1,
            //             prefer N*N+3*N+(N-1)*NB)

            dorgbr('P', N, N, N, VT, LDVT, WORK(ITAUP), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of R in WORK(IU) and computing
            // right singular vectors of R in VT
            // (Workspace: need N*N + BDSPAC)

            dbdsqr(
                'U',
                N,
                N,
                N,
                0,
                S,
                WORK(IE),
                VT,
                LDVT,
                WORK(IU).asMatrix(LDWRKU),
                LDWRKU,
                DUM.asMatrix(1),
                1,
                WORK(IWORK),
                INFO);

            // Multiply Q in A by left singular vectors of R in
            // WORK(IU), storing result in U
            // (Workspace: need N*N)

            dgemm('N', 'N', M, N, N, ONE, A, LDA, WORK(IU).asMatrix(LDWRKU),
                LDWRKU, ZERO, U, LDU);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + N;

            // Compute A=Q*R, copying result to U
            // (Workspace: need 2*N, prefer N + N*NB)

            dgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('L', M, N, A, LDA, U, LDU);

            // Generate Q in U
            // (Workspace: need 2*N, prefer N + N*NB)

            dorgqr(M, N, N, U, LDU, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);

            // Copy R to VT, zeroing out below it

            dlacpy('U', N, N, A, LDA, VT, LDVT);
            if (N > 1) dlaset('L', N - 1, N - 1, ZERO, ZERO, VT(2, 1), LDVT);
            IE = ITAU;
            ITAUQ = IE + N;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Bidiagonalize R in VT
            // (Workspace: need 4*N, prefer 3*N + 2*N*NB)

            dgebrd(N, N, VT, LDVT, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply Q in U by left bidiagonalizing vectors
            // in VT
            // (Workspace: need 3*N + M, prefer 3*N + M*NB)

            dormbr('Q', 'R', 'N', M, N, N, VT, LDVT, WORK(ITAUQ), U, LDU,
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate right bidiagonalizing vectors in VT
            // (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)

            dorgbr('P', N, N, N, VT, LDVT, WORK(ITAUP), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of A in U and computing right
            // singular vectors of A in VT
            // (Workspace: need BDSPAC)

            dbdsqr('U', N, N, M, 0, S, WORK(IE), VT, LDVT, U, LDU,
                DUM.asMatrix(1), 1, WORK(IWORK), INFO);
          }
        }
      } else if (WNTUA) {
        if (WNTVN) {
          // Path 7 (M much larger than N, JOBU='A', JOBVT='N')
          // M left singular vectors to be computed in U and
          // no right singular vectors to be computed

          if (LWORK >= N * N + max(N + M, max(4 * N, BDSPAC))) {
            // Sufficient workspace for a fast algorithm

            IR = 1;
            if (LWORK >= WRKBL + LDA * N) {
              // WORK(IR) is LDA by N

              LDWRKR = LDA;
            } else {
              // WORK(IR) is N by N

              LDWRKR = N;
            }
            ITAU = IR + LDWRKR * N;
            IWORK = ITAU + N;

            // Compute A=Q*R, copying result to U
            // (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)

            dgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('L', M, N, A, LDA, U, LDU);

            // Copy R to WORK(IR), zeroing out below it

            dlacpy('U', N, N, A, LDA, WORK(IR).asMatrix(LDWRKR), LDWRKR);
            dlaset('L', N - 1, N - 1, ZERO, ZERO, WORK(IR + 1).asMatrix(LDWRKR),
                LDWRKR);

            // Generate Q in U
            // (Workspace: need N*N + N + M, prefer N*N + N + M*NB)

            dorgqr(M, M, N, U, LDU, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);
            IE = ITAU;
            ITAUQ = IE + N;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Bidiagonalize R in WORK(IR)
            // (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB)

            dgebrd(N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, S, WORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate left bidiagonalizing vectors in WORK(IR)
            // (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB)

            dorgbr('Q', N, N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUQ),
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of R in WORK(IR)
            // (Workspace: need N*N + BDSPAC)

            dbdsqr(
                'U',
                N,
                0,
                N,
                0,
                S,
                WORK(IE),
                DUM.asMatrix(1),
                1,
                WORK(IR).asMatrix(LDWRKR),
                LDWRKR,
                DUM.asMatrix(1),
                1,
                WORK(IWORK),
                INFO);

            // Multiply Q in U by left singular vectors of R in
            // WORK(IR), storing result in A
            // (Workspace: need N*N)

            dgemm('N', 'N', M, N, N, ONE, U, LDU, WORK(IR).asMatrix(LDWRKR),
                LDWRKR, ZERO, A, LDA);

            // Copy left singular vectors of A from A to U

            dlacpy('F', M, N, A, LDA, U, LDU);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + N;

            // Compute A=Q*R, copying result to U
            // (Workspace: need 2*N, prefer N + N*NB)

            dgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('L', M, N, A, LDA, U, LDU);

            // Generate Q in U
            // (Workspace: need N + M, prefer N + M*NB)

            dorgqr(M, M, N, U, LDU, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);
            IE = ITAU;
            ITAUQ = IE + N;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Zero out below R in A

            if (N > 1) {
              dlaset('L', N - 1, N - 1, ZERO, ZERO, A(2, 1), LDA);
            }

            // Bidiagonalize R in A
            // (Workspace: need 4*N, prefer 3*N + 2*N*NB)

            dgebrd(N, N, A, LDA, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply Q in U by left bidiagonalizing vectors
            // in A
            // (Workspace: need 3*N + M, prefer 3*N + M*NB)

            dormbr('Q', 'R', 'N', M, N, N, A, LDA, WORK(ITAUQ), U, LDU,
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of A in U
            // (Workspace: need BDSPAC)

            dbdsqr('U', N, 0, M, 0, S, WORK(IE), DUM.asMatrix(1), 1, U, LDU,
                DUM.asMatrix(1), 1, WORK(IWORK), INFO);
          }
        } else if (WNTVO) {
          // Path 8 (M much larger than N, JOBU='A', JOBVT='O')
          // M left singular vectors to be computed in U and
          // N right singular vectors to be overwritten on A

          if (LWORK >= 2 * N * N + max(N + M, max(4 * N, BDSPAC))) {
            // Sufficient workspace for a fast algorithm

            IU = 1;
            if (LWORK >= WRKBL + 2 * LDA * N) {
              // WORK(IU) is LDA by N and WORK(IR) is LDA by N

              LDWRKU = LDA;
              IR = IU + LDWRKU * N;
              LDWRKR = LDA;
            } else if (LWORK >= WRKBL + (LDA + N) * N) {
              // WORK(IU) is LDA by N and WORK(IR) is N by N

              LDWRKU = LDA;
              IR = IU + LDWRKU * N;
              LDWRKR = N;
            } else {
              // WORK(IU) is N by N and WORK(IR) is N by N

              LDWRKU = N;
              IR = IU + LDWRKU * N;
              LDWRKR = N;
            }
            ITAU = IR + LDWRKR * N;
            IWORK = ITAU + N;

            // Compute A=Q*R, copying result to U
            // (Workspace: need 2*N*N + 2*N, prefer 2*N*N + N + N*NB)

            dgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('L', M, N, A, LDA, U, LDU);

            // Generate Q in U
            // (Workspace: need 2*N*N + N + M, prefer 2*N*N + N + M*NB)

            dorgqr(M, M, N, U, LDU, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);

            // Copy R to WORK(IU), zeroing out below it

            dlacpy('U', N, N, A, LDA, WORK(IU).asMatrix(LDWRKU), LDWRKU);
            dlaset('L', N - 1, N - 1, ZERO, ZERO, WORK(IU + 1).asMatrix(LDWRKU),
                LDWRKU);
            IE = ITAU;
            ITAUQ = IE + N;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Bidiagonalize R in WORK(IU), copying result to
            // WORK(IR)
            // (Workspace: need 2*N*N + 4*N,
            //             prefer 2*N*N+3*N+2*N*NB)

            dgebrd(N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, S, WORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('U', N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU,
                WORK(IR).asMatrix(LDWRKR), LDWRKR);

            // Generate left bidiagonalizing vectors in WORK(IU)
            // (Workspace: need 2*N*N + 4*N, prefer 2*N*N + 3*N + N*NB)

            dorgbr('Q', N, N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, WORK(ITAUQ),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate right bidiagonalizing vectors in WORK(IR)
            // (Workspace: need 2*N*N + 4*N-1,
            //             prefer 2*N*N+3*N+(N-1)*NB)

            dorgbr('P', N, N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of R in WORK(IU) and computing
            // right singular vectors of R in WORK(IR)
            // (Workspace: need 2*N*N + BDSPAC)

            dbdsqr(
                'U',
                N,
                N,
                N,
                0,
                S,
                WORK(IE),
                WORK(IR).asMatrix(LDWRKR),
                LDWRKR,
                WORK(IU).asMatrix(LDWRKU),
                LDWRKU,
                DUM.asMatrix(1),
                1,
                WORK(IWORK),
                INFO);

            // Multiply Q in U by left singular vectors of R in
            // WORK(IU), storing result in A
            // (Workspace: need N*N)

            dgemm('N', 'N', M, N, N, ONE, U, LDU, WORK(IU).asMatrix(LDWRKU),
                LDWRKU, ZERO, A, LDA);

            // Copy left singular vectors of A from A to U

            dlacpy('F', M, N, A, LDA, U, LDU);

            // Copy right singular vectors of R from WORK(IR) to A

            dlacpy('F', N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, A, LDA);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + N;

            // Compute A=Q*R, copying result to U
            // (Workspace: need 2*N, prefer N + N*NB)

            dgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('L', M, N, A, LDA, U, LDU);

            // Generate Q in U
            // (Workspace: need N + M, prefer N + M*NB)

            dorgqr(M, M, N, U, LDU, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);
            IE = ITAU;
            ITAUQ = IE + N;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Zero out below R in A

            if (N > 1) {
              dlaset('L', N - 1, N - 1, ZERO, ZERO, A(2, 1), LDA);
            }

            // Bidiagonalize R in A
            // (Workspace: need 4*N, prefer 3*N + 2*N*NB)

            dgebrd(N, N, A, LDA, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply Q in U by left bidiagonalizing vectors
            // in A
            // (Workspace: need 3*N + M, prefer 3*N + M*NB)

            dormbr('Q', 'R', 'N', M, N, N, A, LDA, WORK(ITAUQ), U, LDU,
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate right bidiagonalizing vectors in A
            // (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)

            dorgbr('P', N, N, N, A, LDA, WORK(ITAUP), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of A in U and computing right
            // singular vectors of A in A
            // (Workspace: need BDSPAC)

            dbdsqr('U', N, N, M, 0, S, WORK(IE), A, LDA, U, LDU,
                DUM.asMatrix(1), 1, WORK(IWORK), INFO);
          }
        } else if (WNTVAS) {
          // Path 9 (M much larger than N, JOBU='A', JOBVT='S'
          //         or 'A')
          // M left singular vectors to be computed in U and
          // N right singular vectors to be computed in VT

          if (LWORK >= N * N + max(N + M, max(4 * N, BDSPAC))) {
            // Sufficient workspace for a fast algorithm

            IU = 1;
            if (LWORK >= WRKBL + LDA * N) {
              // WORK(IU) is LDA by N

              LDWRKU = LDA;
            } else {
              // WORK(IU) is N by N

              LDWRKU = N;
            }
            ITAU = IU + LDWRKU * N;
            IWORK = ITAU + N;

            // Compute A=Q*R, copying result to U
            // (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)

            dgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('L', M, N, A, LDA, U, LDU);

            // Generate Q in U
            // (Workspace: need N*N + N + M, prefer N*N + N + M*NB)

            dorgqr(M, M, N, U, LDU, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);

            // Copy R to WORK(IU), zeroing out below it

            dlacpy('U', N, N, A, LDA, WORK(IU).asMatrix(LDWRKU), LDWRKU);
            dlaset('L', N - 1, N - 1, ZERO, ZERO, WORK(IU + 1).asMatrix(LDWRKU),
                LDWRKU);
            IE = ITAU;
            ITAUQ = IE + N;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Bidiagonalize R in WORK(IU), copying result to VT
            // (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB)

            dgebrd(N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, S, WORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('U', N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, VT, LDVT);

            // Generate left bidiagonalizing vectors in WORK(IU)
            // (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB)

            dorgbr('Q', N, N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, WORK(ITAUQ),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate right bidiagonalizing vectors in VT
            // (Workspace: need N*N + 4*N-1,
            //             prefer N*N+3*N+(N-1)*NB)

            dorgbr('P', N, N, N, VT, LDVT, WORK(ITAUP), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of R in WORK(IU) and computing
            // right singular vectors of R in VT
            // (Workspace: need N*N + BDSPAC)

            dbdsqr(
                'U',
                N,
                N,
                N,
                0,
                S,
                WORK(IE),
                VT,
                LDVT,
                WORK(IU).asMatrix(LDWRKU),
                LDWRKU,
                DUM.asMatrix(1),
                1,
                WORK(IWORK),
                INFO);

            // Multiply Q in U by left singular vectors of R in
            // WORK(IU), storing result in A
            // (Workspace: need N*N)

            dgemm('N', 'N', M, N, N, ONE, U, LDU, WORK(IU).asMatrix(LDWRKU),
                LDWRKU, ZERO, A, LDA);

            // Copy left singular vectors of A from A to U

            dlacpy('F', M, N, A, LDA, U, LDU);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + N;

            // Compute A=Q*R, copying result to U
            // (Workspace: need 2*N, prefer N + N*NB)

            dgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('L', M, N, A, LDA, U, LDU);

            // Generate Q in U
            // (Workspace: need N + M, prefer N + M*NB)

            dorgqr(M, M, N, U, LDU, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);

            // Copy R from A to VT, zeroing out below it

            dlacpy('U', N, N, A, LDA, VT, LDVT);
            if (N > 1) dlaset('L', N - 1, N - 1, ZERO, ZERO, VT(2, 1), LDVT);
            IE = ITAU;
            ITAUQ = IE + N;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Bidiagonalize R in VT
            // (Workspace: need 4*N, prefer 3*N + 2*N*NB)

            dgebrd(N, N, VT, LDVT, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply Q in U by left bidiagonalizing vectors
            // in VT
            // (Workspace: need 3*N + M, prefer 3*N + M*NB)

            dormbr('Q', 'R', 'N', M, N, N, VT, LDVT, WORK(ITAUQ), U, LDU,
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate right bidiagonalizing vectors in VT
            // (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)

            dorgbr('P', N, N, N, VT, LDVT, WORK(ITAUP), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of A in U and computing right
            // singular vectors of A in VT
            // (Workspace: need BDSPAC)

            dbdsqr('U', N, N, M, 0, S, WORK(IE), VT, LDVT, U, LDU,
                DUM.asMatrix(1), 1, WORK(IWORK), INFO);
          }
        }
      }
    } else {
      // M < MNTHR

      // Path 10 (M at least N, but not much larger)
      // Reduce to bidiagonal form without QR decomposition

      IE = 1;
      ITAUQ = IE + N;
      ITAUP = ITAUQ + N;
      IWORK = ITAUP + N;

      // Bidiagonalize A
      // (Workspace: need 3*N + M, prefer 3*N + (M + N)*NB)

      dgebrd(M, N, A, LDA, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP), WORK(IWORK),
          LWORK - IWORK + 1, IERR);
      if (WNTUAS) {
        // If left singular vectors desired in U, copy result to U
        // and generate left bidiagonalizing vectors in U
        // (Workspace: need 3*N + NCU, prefer 3*N + NCU*NB)

        dlacpy('L', M, N, A, LDA, U, LDU);
        if (WNTUS) NCU = N;
        if (WNTUA) NCU = M;
        dorgbr('Q', M, NCU, N, U, LDU, WORK(ITAUQ), WORK(IWORK),
            LWORK - IWORK + 1, IERR);
      }
      if (WNTVAS) {
        // If right singular vectors desired in VT, copy result to
        // VT and generate right bidiagonalizing vectors in VT
        // (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)

        dlacpy('U', N, N, A, LDA, VT, LDVT);
        dorgbr('P', N, N, N, VT, LDVT, WORK(ITAUP), WORK(IWORK),
            LWORK - IWORK + 1, IERR);
      }
      if (WNTUO) {
        // If left singular vectors desired in A, generate left
        // bidiagonalizing vectors in A
        // (Workspace: need 4*N, prefer 3*N + N*NB)

        dorgbr('Q', M, N, N, A, LDA, WORK(ITAUQ), WORK(IWORK),
            LWORK - IWORK + 1, IERR);
      }
      if (WNTVO) {
        // If right singular vectors desired in A, generate right
        // bidiagonalizing vectors in A
        // (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)

        dorgbr('P', N, N, N, A, LDA, WORK(ITAUP), WORK(IWORK),
            LWORK - IWORK + 1, IERR);
      }
      IWORK = IE + N;
      if (WNTUAS || WNTUO) NRU = M;
      if (WNTUN) NRU = 0;
      if (WNTVAS || WNTVO) NCVT = N;
      if (WNTVN) NCVT = 0;
      if (!WNTUO && !WNTVO) {
        // Perform bidiagonal QR iteration, if desired, computing
        // left singular vectors in U and computing right singular
        // vectors in VT
        // (Workspace: need BDSPAC)

        dbdsqr('U', N, NCVT, NRU, 0, S, WORK(IE), VT, LDVT, U, LDU,
            DUM.asMatrix(1), 1, WORK(IWORK), INFO);
      } else if (!WNTUO && WNTVO) {
        // Perform bidiagonal QR iteration, if desired, computing
        // left singular vectors in U and computing right singular
        // vectors in A
        // (Workspace: need BDSPAC)

        dbdsqr('U', N, NCVT, NRU, 0, S, WORK(IE), A, LDA, U, LDU,
            DUM.asMatrix(1), 1, WORK(IWORK), INFO);
      } else {
        // Perform bidiagonal QR iteration, if desired, computing
        // left singular vectors in A and computing right singular
        // vectors in VT
        // (Workspace: need BDSPAC)

        dbdsqr('U', N, NCVT, NRU, 0, S, WORK(IE), VT, LDVT, A, LDA,
            DUM.asMatrix(1), 1, WORK(IWORK), INFO);
      }
    }
  } else {
    // A has more columns than rows. If A has sufficiently more
    // columns than rows, first reduce using the LQ decomposition (if
    // sufficient workspace available)

    if (N >= MNTHR) {
      if (WNTVN) {
        // Path 1t(N much larger than M, JOBVT='N')
        // No right singular vectors to be computed

        ITAU = 1;
        IWORK = ITAU + M;

        // Compute A=L*Q
        // (Workspace: need 2*M, prefer M + M*NB)

        dgelqf(M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

        // Zero out above L

        dlaset('U', M - 1, M - 1, ZERO, ZERO, A(1, 2), LDA);
        IE = 1;
        ITAUQ = IE + M;
        ITAUP = ITAUQ + M;
        IWORK = ITAUP + M;

        // Bidiagonalize L in A
        // (Workspace: need 4*M, prefer 3*M + 2*M*NB)

        dgebrd(M, M, A, LDA, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP), WORK(IWORK),
            LWORK - IWORK + 1, IERR);
        if (WNTUO || WNTUAS) {
          // If left singular vectors desired, generate Q
          // (Workspace: need 4*M, prefer 3*M + M*NB)

          dorgbr('Q', M, M, M, A, LDA, WORK(ITAUQ), WORK(IWORK),
              LWORK - IWORK + 1, IERR);
        }
        IWORK = IE + M;
        NRU = 0;
        if (WNTUO || WNTUAS) NRU = M;

        // Perform bidiagonal QR iteration, computing left singular
        // vectors of A in A if desired
        // (Workspace: need BDSPAC)

        dbdsqr('U', M, 0, NRU, 0, S, WORK(IE), DUM.asMatrix(1), 1, A, LDA,
            DUM.asMatrix(1), 1, WORK(IWORK), INFO);

        // If left singular vectors desired in U, copy them there

        if (WNTUAS) dlacpy('F', M, M, A, LDA, U, LDU);
      } else if (WNTVO && WNTUN) {
        // Path 2t(N much larger than M, JOBU='N', JOBVT='O')
        // M right singular vectors to be overwritten on A and
        // no left singular vectors to be computed

        if (LWORK >= M * M + max(4 * M, BDSPAC)) {
          // Sufficient workspace for a fast algorithm

          IR = 1;
          if (LWORK >= max(WRKBL, LDA * N + M) + LDA * M) {
            // WORK(IU) is LDA by N and WORK(IR) is LDA by M

            LDWRKU = LDA;
            CHUNK = N;
            LDWRKR = LDA;
          } else if (LWORK >= max(WRKBL, LDA * N + M) + M * M) {
            // WORK(IU) is LDA by N and WORK(IR) is M by M

            LDWRKU = LDA;
            CHUNK = N;
            LDWRKR = M;
          } else {
            // WORK(IU) is M by CHUNK and WORK(IR) is M by M

            LDWRKU = M;
            CHUNK = (LWORK - M * M - M) ~/ M;
            LDWRKR = M;
          }
          ITAU = IR + LDWRKR * M;
          IWORK = ITAU + M;

          // Compute A=L*Q
          // (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)

          dgelqf(
              M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Copy L to WORK(IR) and zero out above it

          dlacpy('L', M, M, A, LDA, WORK(IR).asMatrix(LDWRKR), LDWRKR);
          dlaset('U', M - 1, M - 1, ZERO, ZERO,
              WORK(IR + LDWRKR).asMatrix(LDWRKR), LDWRKR);

          // Generate Q in A
          // (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)

          dorglq(M, N, M, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
              IERR);
          IE = ITAU;
          ITAUQ = IE + M;
          ITAUP = ITAUQ + M;
          IWORK = ITAUP + M;

          // Bidiagonalize L in WORK(IR)
          // (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB)

          dgebrd(M, M, WORK(IR).asMatrix(LDWRKR), LDWRKR, S, WORK(IE),
              WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Generate right vectors bidiagonalizing L
          // (Workspace: need M*M + 4*M-1, prefer M*M + 3*M + (M-1)*NB)

          dorgbr('P', M, M, M, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUP),
              WORK(IWORK), LWORK - IWORK + 1, IERR);
          IWORK = IE + M;

          // Perform bidiagonal QR iteration, computing right
          // singular vectors of L in WORK(IR)
          // (Workspace: need M*M + BDSPAC)

          dbdsqr(
              'U',
              M,
              M,
              0,
              0,
              S,
              WORK(IE),
              WORK(IR).asMatrix(LDWRKR),
              LDWRKR,
              DUM.asMatrix(1),
              1,
              DUM.asMatrix(1),
              1,
              WORK(IWORK),
              INFO);
          IU = IE + M;

          // Multiply right singular vectors of L in WORK(IR) by Q
          // in A, storing result in WORK(IU) and copying to A
          // (Workspace: need M*M + 2*M, prefer M*M + M*N + M)

          for (I = 1; I <= N; I += CHUNK) {
            BLK = min(N - I + 1, CHUNK);
            dgemm('N', 'N', M, BLK, M, ONE, WORK(IR).asMatrix(LDWRKR), LDWRKR,
                A(1, I), LDA, ZERO, WORK(IU).asMatrix(LDWRKU), LDWRKU);
            dlacpy(
                'F', M, BLK, WORK(IU).asMatrix(LDWRKU), LDWRKU, A(1, I), LDA);
          }
        } else {
          // Insufficient workspace for a fast algorithm

          IE = 1;
          ITAUQ = IE + M;
          ITAUP = ITAUQ + M;
          IWORK = ITAUP + M;

          // Bidiagonalize A
          // (Workspace: need 3*M + N, prefer 3*M + (M + N)*NB)

          dgebrd(M, N, A, LDA, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP),
              WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Generate right vectors bidiagonalizing A
          // (Workspace: need 4*M, prefer 3*M + M*NB)

          dorgbr('P', M, N, M, A, LDA, WORK(ITAUP), WORK(IWORK),
              LWORK - IWORK + 1, IERR);
          IWORK = IE + M;

          // Perform bidiagonal QR iteration, computing right
          // singular vectors of A in A
          // (Workspace: need BDSPAC)

          dbdsqr('L', M, N, 0, 0, S, WORK(IE), A, LDA, DUM.asMatrix(1), 1,
              DUM.asMatrix(1), 1, WORK(IWORK), INFO);
        }
      } else if (WNTVO && WNTUAS) {
        // Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O')
        // M right singular vectors to be overwritten on A and
        // M left singular vectors to be computed in U

        if (LWORK >= M * M + max(4 * M, BDSPAC)) {
          // Sufficient workspace for a fast algorithm

          IR = 1;
          if (LWORK >= max(WRKBL, LDA * N + M) + LDA * M) {
            // WORK(IU) is LDA by N and WORK(IR) is LDA by M

            LDWRKU = LDA;
            CHUNK = N;
            LDWRKR = LDA;
          } else if (LWORK >= max(WRKBL, LDA * N + M) + M * M) {
            // WORK(IU) is LDA by N and WORK(IR) is M by M

            LDWRKU = LDA;
            CHUNK = N;
            LDWRKR = M;
          } else {
            // WORK(IU) is M by CHUNK and WORK(IR) is M by M

            LDWRKU = M;
            CHUNK = (LWORK - M * M - M) ~/ M;
            LDWRKR = M;
          }
          ITAU = IR + LDWRKR * M;
          IWORK = ITAU + M;

          // Compute A=L*Q
          // (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)

          dgelqf(
              M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Copy L to U, zeroing about above it

          dlacpy('L', M, M, A, LDA, U, LDU);
          dlaset('U', M - 1, M - 1, ZERO, ZERO, U(1, 2), LDU);

          // Generate Q in A
          // (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)

          dorglq(M, N, M, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
              IERR);
          IE = ITAU;
          ITAUQ = IE + M;
          ITAUP = ITAUQ + M;
          IWORK = ITAUP + M;

          // Bidiagonalize L in U, copying result to WORK(IR)
          // (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB)

          dgebrd(M, M, U, LDU, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP),
              WORK(IWORK), LWORK - IWORK + 1, IERR);
          dlacpy('U', M, M, U, LDU, WORK(IR).asMatrix(LDWRKR), LDWRKR);

          // Generate right vectors bidiagonalizing L in WORK(IR)
          // (Workspace: need M*M + 4*M-1, prefer M*M + 3*M + (M-1)*NB)

          dorgbr('P', M, M, M, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUP),
              WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Generate left vectors bidiagonalizing L in U
          // (Workspace: need M*M + 4*M, prefer M*M + 3*M + M*NB)

          dorgbr('Q', M, M, M, U, LDU, WORK(ITAUQ), WORK(IWORK),
              LWORK - IWORK + 1, IERR);
          IWORK = IE + M;

          // Perform bidiagonal QR iteration, computing left
          // singular vectors of L in U, and computing right
          // singular vectors of L in WORK(IR)
          // (Workspace: need M*M + BDSPAC)

          dbdsqr('U', M, M, M, 0, S, WORK(IE), WORK(IR).asMatrix(LDWRKR),
              LDWRKR, U, LDU, DUM.asMatrix(1), 1, WORK(IWORK), INFO);
          IU = IE + M;

          // Multiply right singular vectors of L in WORK(IR) by Q
          // in A, storing result in WORK(IU) and copying to A
          // (Workspace: need M*M + 2*M, prefer M*M + M*N + M))

          for (I = 1; I <= N; I += CHUNK) {
            BLK = min(N - I + 1, CHUNK);
            dgemm('N', 'N', M, BLK, M, ONE, WORK(IR).asMatrix(LDWRKR), LDWRKR,
                A(1, I), LDA, ZERO, WORK(IU).asMatrix(LDWRKU), LDWRKU);
            dlacpy(
                'F', M, BLK, WORK(IU).asMatrix(LDWRKU), LDWRKU, A(1, I), LDA);
          }
        } else {
          // Insufficient workspace for a fast algorithm

          ITAU = 1;
          IWORK = ITAU + M;

          // Compute A=L*Q
          // (Workspace: need 2*M, prefer M + M*NB)

          dgelqf(
              M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Copy L to U, zeroing out above it

          dlacpy('L', M, M, A, LDA, U, LDU);
          dlaset('U', M - 1, M - 1, ZERO, ZERO, U(1, 2), LDU);

          // Generate Q in A
          // (Workspace: need 2*M, prefer M + M*NB)

          dorglq(M, N, M, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
              IERR);
          IE = ITAU;
          ITAUQ = IE + M;
          ITAUP = ITAUQ + M;
          IWORK = ITAUP + M;

          // Bidiagonalize L in U
          // (Workspace: need 4*M, prefer 3*M + 2*M*NB)

          dgebrd(M, M, U, LDU, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP),
              WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Multiply right vectors bidiagonalizing L by Q in A
          // (Workspace: need 3*M + N, prefer 3*M + N*NB)

          dormbr('P', 'L', 'T', M, N, M, U, LDU, WORK(ITAUP), A, LDA,
              WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Generate left vectors bidiagonalizing L in U
          // (Workspace: need 4*M, prefer 3*M + M*NB)

          dorgbr('Q', M, M, M, U, LDU, WORK(ITAUQ), WORK(IWORK),
              LWORK - IWORK + 1, IERR);
          IWORK = IE + M;

          // Perform bidiagonal QR iteration, computing left
          // singular vectors of A in U and computing right
          // singular vectors of A in A
          // (Workspace: need BDSPAC)

          dbdsqr('U', M, N, M, 0, S, WORK(IE), A, LDA, U, LDU, DUM.asMatrix(1),
              1, WORK(IWORK), INFO);
        }
      } else if (WNTVS) {
        if (WNTUN) {
          // Path 4t(N much larger than M, JOBU='N', JOBVT='S')
          // M right singular vectors to be computed in VT and
          // no left singular vectors to be computed

          if (LWORK >= M * M + max(4 * M, BDSPAC)) {
            // Sufficient workspace for a fast algorithm

            IR = 1;
            if (LWORK >= WRKBL + LDA * M) {
              // WORK(IR) is LDA by M

              LDWRKR = LDA;
            } else {
              // WORK(IR) is M by M

              LDWRKR = M;
            }
            ITAU = IR + LDWRKR * M;
            IWORK = ITAU + M;

            // Compute A=L*Q
            // (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)

            dgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Copy L to WORK(IR), zeroing out above it

            dlacpy('L', M, M, A, LDA, WORK(IR).asMatrix(LDWRKR), LDWRKR);
            dlaset('U', M - 1, M - 1, ZERO, ZERO,
                WORK(IR + LDWRKR).asMatrix(LDWRKR), LDWRKR);

            // Generate Q in A
            // (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)

            dorglq(M, N, M, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);
            IE = ITAU;
            ITAUQ = IE + M;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Bidiagonalize L in WORK(IR)
            // (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB)

            dgebrd(M, M, WORK(IR).asMatrix(LDWRKR), LDWRKR, S, WORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate right vectors bidiagonalizing L in
            // WORK(IR)
            // (Workspace: need M*M + 4*M, prefer M*M + 3*M + (M-1)*NB)

            dorgbr('P', M, M, M, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IWORK = IE + M;

            // Perform bidiagonal QR iteration, computing right
            // singular vectors of L in WORK(IR)
            // (Workspace: need M*M + BDSPAC)

            dbdsqr(
                'U',
                M,
                M,
                0,
                0,
                S,
                WORK(IE),
                WORK(IR).asMatrix(LDWRKR),
                LDWRKR,
                DUM.asMatrix(1),
                1,
                DUM.asMatrix(1),
                1,
                WORK(IWORK),
                INFO);

            // Multiply right singular vectors of L in WORK(IR) by
            // Q in A, storing result in VT
            // (Workspace: need M*M)

            dgemm('N', 'N', M, N, M, ONE, WORK(IR).asMatrix(LDWRKR), LDWRKR, A,
                LDA, ZERO, VT, LDVT);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + M;

            // Compute A=L*Q
            // (Workspace: need 2*M, prefer M + M*NB)

            dgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Copy result to VT

            dlacpy('U', M, N, A, LDA, VT, LDVT);

            // Generate Q in VT
            // (Workspace: need 2*M, prefer M + M*NB)

            dorglq(M, N, M, VT, LDVT, WORK(ITAU), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IE = ITAU;
            ITAUQ = IE + M;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Zero out above L in A

            dlaset('U', M - 1, M - 1, ZERO, ZERO, A(1, 2), LDA);

            // Bidiagonalize L in A
            // (Workspace: need 4*M, prefer 3*M + 2*M*NB)

            dgebrd(M, M, A, LDA, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply right vectors bidiagonalizing L by Q in VT
            // (Workspace: need 3*M + N, prefer 3*M + N*NB)

            dormbr('P', 'L', 'T', M, N, M, A, LDA, WORK(ITAUP), VT, LDVT,
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IWORK = IE + M;

            // Perform bidiagonal QR iteration, computing right
            // singular vectors of A in VT
            // (Workspace: need BDSPAC)

            dbdsqr('U', M, N, 0, 0, S, WORK(IE), VT, LDVT, DUM.asMatrix(1), 1,
                DUM.asMatrix(1), 1, WORK(IWORK), INFO);
          }
        } else if (WNTUO) {
          // Path 5t(N much larger than M, JOBU='O', JOBVT='S')
          // M right singular vectors to be computed in VT and
          // M left singular vectors to be overwritten on A

          if (LWORK >= 2 * M * M + max(4 * M, BDSPAC)) {
            // Sufficient workspace for a fast algorithm

            IU = 1;
            if (LWORK >= WRKBL + 2 * LDA * M) {
              // WORK(IU) is LDA by M and WORK(IR) is LDA by M

              LDWRKU = LDA;
              IR = IU + LDWRKU * M;
              LDWRKR = LDA;
            } else if (LWORK >= WRKBL + (LDA + M) * M) {
              // WORK(IU) is LDA by M and WORK(IR) is M by M

              LDWRKU = LDA;
              IR = IU + LDWRKU * M;
              LDWRKR = M;
            } else {
              // WORK(IU) is M by M and WORK(IR) is M by M

              LDWRKU = M;
              IR = IU + LDWRKU * M;
              LDWRKR = M;
            }
            ITAU = IR + LDWRKR * M;
            IWORK = ITAU + M;

            // Compute A=L*Q
            // (Workspace: need 2*M*M + 2*M, prefer 2*M*M + M + M*NB)

            dgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Copy L to WORK(IU), zeroing out below it

            dlacpy('L', M, M, A, LDA, WORK(IU).asMatrix(LDWRKU), LDWRKU);
            dlaset('U', M - 1, M - 1, ZERO, ZERO,
                WORK(IU + LDWRKU).asMatrix(LDWRKU), LDWRKU);

            // Generate Q in A
            // (Workspace: need 2*M*M + 2*M, prefer 2*M*M + M + M*NB)

            dorglq(M, N, M, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);
            IE = ITAU;
            ITAUQ = IE + M;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Bidiagonalize L in WORK(IU), copying result to
            // WORK(IR)
            // (Workspace: need 2*M*M + 4*M,
            //             prefer 2*M*M+3*M+2*M*NB)

            dgebrd(M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU, S, WORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('L', M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU,
                WORK(IR).asMatrix(LDWRKR), LDWRKR);

            // Generate right bidiagonalizing vectors in WORK(IU)
            // (Workspace: need 2*M*M + 4*M-1,
            //             prefer 2*M*M+3*M+(M-1)*NB)

            dorgbr('P', M, M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU, WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate left bidiagonalizing vectors in WORK(IR)
            // (Workspace: need 2*M*M + 4*M, prefer 2*M*M + 3*M + M*NB)

            dorgbr('Q', M, M, M, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUQ),
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IWORK = IE + M;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of L in WORK(IR) and computing
            // right singular vectors of L in WORK(IU)
            // (Workspace: need 2*M*M + BDSPAC)

            dbdsqr(
                'U',
                M,
                M,
                M,
                0,
                S,
                WORK(IE),
                WORK(IU).asMatrix(LDWRKU),
                LDWRKU,
                WORK(IR).asMatrix(LDWRKR),
                LDWRKR,
                DUM.asMatrix(1),
                1,
                WORK(IWORK),
                INFO);

            // Multiply right singular vectors of L in WORK(IU) by
            // Q in A, storing result in VT
            // (Workspace: need M*M)

            dgemm('N', 'N', M, N, M, ONE, WORK(IU).asMatrix(LDWRKU), LDWRKU, A,
                LDA, ZERO, VT, LDVT);

            // Copy left singular vectors of L to A
            // (Workspace: need M*M)

            dlacpy('F', M, M, WORK(IR).asMatrix(LDWRKR), LDWRKR, A, LDA);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + M;

            // Compute A=L*Q, copying result to VT
            // (Workspace: need 2*M, prefer M + M*NB)

            dgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('U', M, N, A, LDA, VT, LDVT);

            // Generate Q in VT
            // (Workspace: need 2*M, prefer M + M*NB)

            dorglq(M, N, M, VT, LDVT, WORK(ITAU), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IE = ITAU;
            ITAUQ = IE + M;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Zero out above L in A

            dlaset('U', M - 1, M - 1, ZERO, ZERO, A(1, 2), LDA);

            // Bidiagonalize L in A
            // (Workspace: need 4*M, prefer 3*M + 2*M*NB)

            dgebrd(M, M, A, LDA, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply right vectors bidiagonalizing L by Q in VT
            // (Workspace: need 3*M + N, prefer 3*M + N*NB)

            dormbr('P', 'L', 'T', M, N, M, A, LDA, WORK(ITAUP), VT, LDVT,
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate left bidiagonalizing vectors of L in A
            // (Workspace: need 4*M, prefer 3*M + M*NB)

            dorgbr('Q', M, M, M, A, LDA, WORK(ITAUQ), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IWORK = IE + M;

            // Perform bidiagonal QR iteration, compute left
            // singular vectors of A in A and compute right
            // singular vectors of A in VT
            // (Workspace: need BDSPAC)

            dbdsqr('U', M, N, M, 0, S, WORK(IE), VT, LDVT, A, LDA,
                DUM.asMatrix(1), 1, WORK(IWORK), INFO);
          }
        } else if (WNTUAS) {
          // Path 6t(N much larger than M, JOBU='S' or 'A',
          //         JOBVT='S')
          // M right singular vectors to be computed in VT and
          // M left singular vectors to be computed in U

          if (LWORK >= M * M + max(4 * M, BDSPAC)) {
            // Sufficient workspace for a fast algorithm

            IU = 1;
            if (LWORK >= WRKBL + LDA * M) {
              // WORK(IU) is LDA by N

              LDWRKU = LDA;
            } else {
              // WORK(IU) is LDA by M

              LDWRKU = M;
            }
            ITAU = IU + LDWRKU * M;
            IWORK = ITAU + M;

            // Compute A=L*Q
            // (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)

            dgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Copy L to WORK(IU), zeroing out above it

            dlacpy('L', M, M, A, LDA, WORK(IU).asMatrix(LDWRKU), LDWRKU);
            dlaset('U', M - 1, M - 1, ZERO, ZERO,
                WORK(IU + LDWRKU).asMatrix(LDWRKU), LDWRKU);

            // Generate Q in A
            // (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)

            dorglq(M, N, M, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);
            IE = ITAU;
            ITAUQ = IE + M;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Bidiagonalize L in WORK(IU), copying result to U
            // (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB)

            dgebrd(M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU, S, WORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('L', M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU, U, LDU);

            // Generate right bidiagonalizing vectors in WORK(IU)
            // (Workspace: need M*M + 4*M-1,
            //             prefer M*M+3*M+(M-1)*NB)

            dorgbr('P', M, M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU, WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate left bidiagonalizing vectors in U
            // (Workspace: need M*M + 4*M, prefer M*M + 3*M + M*NB)

            dorgbr('Q', M, M, M, U, LDU, WORK(ITAUQ), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IWORK = IE + M;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of L in U and computing right
            // singular vectors of L in WORK(IU)
            // (Workspace: need M*M + BDSPAC)

            dbdsqr('U', M, M, M, 0, S, WORK(IE), WORK(IU).asMatrix(LDWRKU),
                LDWRKU, U, LDU, DUM.asMatrix(1), 1, WORK(IWORK), INFO);

            // Multiply right singular vectors of L in WORK(IU) by
            // Q in A, storing result in VT
            // (Workspace: need M*M)

            dgemm('N', 'N', M, N, M, ONE, WORK(IU).asMatrix(LDWRKU), LDWRKU, A,
                LDA, ZERO, VT, LDVT);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + M;

            // Compute A=L*Q, copying result to VT
            // (Workspace: need 2*M, prefer M + M*NB)

            dgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('U', M, N, A, LDA, VT, LDVT);

            // Generate Q in VT
            // (Workspace: need 2*M, prefer M + M*NB)

            dorglq(M, N, M, VT, LDVT, WORK(ITAU), WORK(IWORK),
                LWORK - IWORK + 1, IERR);

            // Copy L to U, zeroing out above it

            dlacpy('L', M, M, A, LDA, U, LDU);
            dlaset('U', M - 1, M - 1, ZERO, ZERO, U(1, 2), LDU);
            IE = ITAU;
            ITAUQ = IE + M;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Bidiagonalize L in U
            // (Workspace: need 4*M, prefer 3*M + 2*M*NB)

            dgebrd(M, M, U, LDU, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply right bidiagonalizing vectors in U by Q
            // in VT
            // (Workspace: need 3*M + N, prefer 3*M + N*NB)

            dormbr('P', 'L', 'T', M, N, M, U, LDU, WORK(ITAUP), VT, LDVT,
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate left bidiagonalizing vectors in U
            // (Workspace: need 4*M, prefer 3*M + M*NB)

            dorgbr('Q', M, M, M, U, LDU, WORK(ITAUQ), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IWORK = IE + M;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of A in U and computing right
            // singular vectors of A in VT
            // (Workspace: need BDSPAC)

            dbdsqr('U', M, N, M, 0, S, WORK(IE), VT, LDVT, U, LDU,
                DUM.asMatrix(1), 1, WORK(IWORK), INFO);
          }
        }
      } else if (WNTVA) {
        if (WNTUN) {
          // Path 7t(N much larger than M, JOBU='N', JOBVT='A')
          // N right singular vectors to be computed in VT and
          // no left singular vectors to be computed

          if (LWORK >= M * M + max(N + M, max(4 * M, BDSPAC))) {
            // Sufficient workspace for a fast algorithm

            IR = 1;
            if (LWORK >= WRKBL + LDA * M) {
              // WORK(IR) is LDA by M

              LDWRKR = LDA;
            } else {
              // WORK(IR) is M by M

              LDWRKR = M;
            }
            ITAU = IR + LDWRKR * M;
            IWORK = ITAU + M;

            // Compute A=L*Q, copying result to VT
            // (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)

            dgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('U', M, N, A, LDA, VT, LDVT);

            // Copy L to WORK(IR), zeroing out above it

            dlacpy('L', M, M, A, LDA, WORK(IR).asMatrix(LDWRKR), LDWRKR);
            dlaset('U', M - 1, M - 1, ZERO, ZERO,
                WORK(IR + LDWRKR).asMatrix(LDWRKR), LDWRKR);

            // Generate Q in VT
            // (Workspace: need M*M + M + N, prefer M*M + M + N*NB)

            dorglq(N, N, M, VT, LDVT, WORK(ITAU), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IE = ITAU;
            ITAUQ = IE + M;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Bidiagonalize L in WORK(IR)
            // (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB)

            dgebrd(M, M, WORK(IR).asMatrix(LDWRKR), LDWRKR, S, WORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate right bidiagonalizing vectors in WORK(IR)
            // (Workspace: need M*M + 4*M-1,
            //             prefer M*M+3*M+(M-1)*NB)

            dorgbr('P', M, M, M, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IWORK = IE + M;

            // Perform bidiagonal QR iteration, computing right
            // singular vectors of L in WORK(IR)
            // (Workspace: need M*M + BDSPAC)

            dbdsqr(
                'U',
                M,
                M,
                0,
                0,
                S,
                WORK(IE),
                WORK(IR).asMatrix(LDWRKR),
                LDWRKR,
                DUM.asMatrix(1),
                1,
                DUM.asMatrix(1),
                1,
                WORK(IWORK),
                INFO);

            // Multiply right singular vectors of L in WORK(IR) by
            // Q in VT, storing result in A
            // (Workspace: need M*M)

            dgemm('N', 'N', M, N, M, ONE, WORK(IR).asMatrix(LDWRKR), LDWRKR, VT,
                LDVT, ZERO, A, LDA);

            // Copy right singular vectors of A from A to VT

            dlacpy('F', M, N, A, LDA, VT, LDVT);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + M;

            // Compute A=L*Q, copying result to VT
            // (Workspace: need 2*M, prefer M + M*NB)

            dgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('U', M, N, A, LDA, VT, LDVT);

            // Generate Q in VT
            // (Workspace: need M + N, prefer M + N*NB)

            dorglq(N, N, M, VT, LDVT, WORK(ITAU), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IE = ITAU;
            ITAUQ = IE + M;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Zero out above L in A

            dlaset('U', M - 1, M - 1, ZERO, ZERO, A(1, 2), LDA);

            // Bidiagonalize L in A
            // (Workspace: need 4*M, prefer 3*M + 2*M*NB)

            dgebrd(M, M, A, LDA, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply right bidiagonalizing vectors in A by Q
            // in VT
            // (Workspace: need 3*M + N, prefer 3*M + N*NB)

            dormbr('P', 'L', 'T', M, N, M, A, LDA, WORK(ITAUP), VT, LDVT,
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IWORK = IE + M;

            // Perform bidiagonal QR iteration, computing right
            // singular vectors of A in VT
            // (Workspace: need BDSPAC)

            dbdsqr('U', M, N, 0, 0, S, WORK(IE), VT, LDVT, DUM.asMatrix(1), 1,
                DUM.asMatrix(1), 1, WORK(IWORK), INFO);
          }
        } else if (WNTUO) {
          // Path 8t(N much larger than M, JOBU='O', JOBVT='A')
          // N right singular vectors to be computed in VT and
          // M left singular vectors to be overwritten on A

          if (LWORK >= 2 * M * M + max(N + M, max(4 * M, BDSPAC))) {
            // Sufficient workspace for a fast algorithm

            IU = 1;
            if (LWORK >= WRKBL + 2 * LDA * M) {
              // WORK(IU) is LDA by M and WORK(IR) is LDA by M

              LDWRKU = LDA;
              IR = IU + LDWRKU * M;
              LDWRKR = LDA;
            } else if (LWORK >= WRKBL + (LDA + M) * M) {
              // WORK(IU) is LDA by M and WORK(IR) is M by M

              LDWRKU = LDA;
              IR = IU + LDWRKU * M;
              LDWRKR = M;
            } else {
              // WORK(IU) is M by M and WORK(IR) is M by M

              LDWRKU = M;
              IR = IU + LDWRKU * M;
              LDWRKR = M;
            }
            ITAU = IR + LDWRKR * M;
            IWORK = ITAU + M;

            // Compute A=L*Q, copying result to VT
            // (Workspace: need 2*M*M + 2*M, prefer 2*M*M + M + M*NB)

            dgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('U', M, N, A, LDA, VT, LDVT);

            // Generate Q in VT
            // (Workspace: need 2*M*M + M + N, prefer 2*M*M + M + N*NB)

            dorglq(N, N, M, VT, LDVT, WORK(ITAU), WORK(IWORK),
                LWORK - IWORK + 1, IERR);

            // Copy L to WORK(IU), zeroing out above it

            dlacpy('L', M, M, A, LDA, WORK(IU).asMatrix(LDWRKU), LDWRKU);
            dlaset('U', M - 1, M - 1, ZERO, ZERO,
                WORK(IU + LDWRKU).asMatrix(LDWRKU), LDWRKU);
            IE = ITAU;
            ITAUQ = IE + M;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Bidiagonalize L in WORK(IU), copying result to
            // WORK(IR)
            // (Workspace: need 2*M*M + 4*M,
            //             prefer 2*M*M+3*M+2*M*NB)

            dgebrd(M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU, S, WORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('L', M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU,
                WORK(IR).asMatrix(LDWRKR), LDWRKR);

            // Generate right bidiagonalizing vectors in WORK(IU)
            // (Workspace: need 2*M*M + 4*M-1,
            //             prefer 2*M*M+3*M+(M-1)*NB)

            dorgbr('P', M, M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU, WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate left bidiagonalizing vectors in WORK(IR)
            // (Workspace: need 2*M*M + 4*M, prefer 2*M*M + 3*M + M*NB)

            dorgbr('Q', M, M, M, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUQ),
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IWORK = IE + M;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of L in WORK(IR) and computing
            // right singular vectors of L in WORK(IU)
            // (Workspace: need 2*M*M + BDSPAC)

            dbdsqr(
                'U',
                M,
                M,
                M,
                0,
                S,
                WORK(IE),
                WORK(IU).asMatrix(LDWRKU),
                LDWRKU,
                WORK(IR).asMatrix(LDWRKR),
                LDWRKR,
                DUM.asMatrix(1),
                1,
                WORK(IWORK),
                INFO);

            // Multiply right singular vectors of L in WORK(IU) by
            // Q in VT, storing result in A
            // (Workspace: need M*M)

            dgemm('N', 'N', M, N, M, ONE, WORK(IU).asMatrix(LDWRKU), LDWRKU, VT,
                LDVT, ZERO, A, LDA);

            // Copy right singular vectors of A from A to VT

            dlacpy('F', M, N, A, LDA, VT, LDVT);

            // Copy left singular vectors of A from WORK(IR) to A

            dlacpy('F', M, M, WORK(IR).asMatrix(LDWRKR), LDWRKR, A, LDA);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + M;

            // Compute A=L*Q, copying result to VT
            // (Workspace: need 2*M, prefer M + M*NB)

            dgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('U', M, N, A, LDA, VT, LDVT);

            // Generate Q in VT
            // (Workspace: need M + N, prefer M + N*NB)

            dorglq(N, N, M, VT, LDVT, WORK(ITAU), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IE = ITAU;
            ITAUQ = IE + M;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Zero out above L in A

            dlaset('U', M - 1, M - 1, ZERO, ZERO, A(1, 2), LDA);

            // Bidiagonalize L in A
            // (Workspace: need 4*M, prefer 3*M + 2*M*NB)

            dgebrd(M, M, A, LDA, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply right bidiagonalizing vectors in A by Q
            // in VT
            // (Workspace: need 3*M + N, prefer 3*M + N*NB)

            dormbr('P', 'L', 'T', M, N, M, A, LDA, WORK(ITAUP), VT, LDVT,
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate left bidiagonalizing vectors in A
            // (Workspace: need 4*M, prefer 3*M + M*NB)

            dorgbr('Q', M, M, M, A, LDA, WORK(ITAUQ), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IWORK = IE + M;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of A in A and computing right
            // singular vectors of A in VT
            // (Workspace: need BDSPAC)

            dbdsqr('U', M, N, M, 0, S, WORK(IE), VT, LDVT, A, LDA,
                DUM.asMatrix(1), 1, WORK(IWORK), INFO);
          }
        } else if (WNTUAS) {
          // Path 9t(N much larger than M, JOBU='S' or 'A',
          //         JOBVT='A')
          // N right singular vectors to be computed in VT and
          // M left singular vectors to be computed in U

          if (LWORK >= M * M + max(N + M, max(4 * M, BDSPAC))) {
            // Sufficient workspace for a fast algorithm

            IU = 1;
            if (LWORK >= WRKBL + LDA * M) {
              // WORK(IU) is LDA by M

              LDWRKU = LDA;
            } else {
              // WORK(IU) is M by M

              LDWRKU = M;
            }
            ITAU = IU + LDWRKU * M;
            IWORK = ITAU + M;

            // Compute A=L*Q, copying result to VT
            // (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)

            dgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('U', M, N, A, LDA, VT, LDVT);

            // Generate Q in VT
            // (Workspace: need M*M + M + N, prefer M*M + M + N*NB)

            dorglq(N, N, M, VT, LDVT, WORK(ITAU), WORK(IWORK),
                LWORK - IWORK + 1, IERR);

            // Copy L to WORK(IU), zeroing out above it

            dlacpy('L', M, M, A, LDA, WORK(IU).asMatrix(LDWRKU), LDWRKU);
            dlaset('U', M - 1, M - 1, ZERO, ZERO,
                WORK(IU + LDWRKU).asMatrix(LDWRKU), LDWRKU);
            IE = ITAU;
            ITAUQ = IE + M;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Bidiagonalize L in WORK(IU), copying result to U
            // (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB)

            dgebrd(M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU, S, WORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('L', M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU, U, LDU);

            // Generate right bidiagonalizing vectors in WORK(IU)
            // (Workspace: need M*M + 4*M, prefer M*M + 3*M + (M-1)*NB)

            dorgbr('P', M, M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU, WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate left bidiagonalizing vectors in U
            // (Workspace: need M*M + 4*M, prefer M*M + 3*M + M*NB)

            dorgbr('Q', M, M, M, U, LDU, WORK(ITAUQ), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IWORK = IE + M;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of L in U and computing right
            // singular vectors of L in WORK(IU)
            // (Workspace: need M*M + BDSPAC)

            dbdsqr('U', M, M, M, 0, S, WORK(IE), WORK(IU).asMatrix(LDWRKU),
                LDWRKU, U, LDU, DUM.asMatrix(1), 1, WORK(IWORK), INFO);

            // Multiply right singular vectors of L in WORK(IU) by
            // Q in VT, storing result in A
            // (Workspace: need M*M)

            dgemm('N', 'N', M, N, M, ONE, WORK(IU).asMatrix(LDWRKU), LDWRKU, VT,
                LDVT, ZERO, A, LDA);

            // Copy right singular vectors of A from A to VT

            dlacpy('F', M, N, A, LDA, VT, LDVT);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + M;

            // Compute A=L*Q, copying result to VT
            // (Workspace: need 2*M, prefer M + M*NB)

            dgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            dlacpy('U', M, N, A, LDA, VT, LDVT);

            // Generate Q in VT
            // (Workspace: need M + N, prefer M + N*NB)

            dorglq(N, N, M, VT, LDVT, WORK(ITAU), WORK(IWORK),
                LWORK - IWORK + 1, IERR);

            // Copy L to U, zeroing out above it

            dlacpy('L', M, M, A, LDA, U, LDU);
            dlaset('U', M - 1, M - 1, ZERO, ZERO, U(1, 2), LDU);
            IE = ITAU;
            ITAUQ = IE + M;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Bidiagonalize L in U
            // (Workspace: need 4*M, prefer 3*M + 2*M*NB)

            dgebrd(M, M, U, LDU, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply right bidiagonalizing vectors in U by Q
            // in VT
            // (Workspace: need 3*M + N, prefer 3*M + N*NB)

            dormbr('P', 'L', 'T', M, N, M, U, LDU, WORK(ITAUP), VT, LDVT,
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate left bidiagonalizing vectors in U
            // (Workspace: need 4*M, prefer 3*M + M*NB)

            dorgbr('Q', M, M, M, U, LDU, WORK(ITAUQ), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IWORK = IE + M;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of A in U and computing right
            // singular vectors of A in VT
            // (Workspace: need BDSPAC)

            dbdsqr('U', M, N, M, 0, S, WORK(IE), VT, LDVT, U, LDU,
                DUM.asMatrix(1), 1, WORK(IWORK), INFO);
          }
        }
      }
    } else {
      // N < MNTHR

      // Path 10t(N greater than M, but not much larger)
      // Reduce to bidiagonal form without LQ decomposition

      IE = 1;
      ITAUQ = IE + M;
      ITAUP = ITAUQ + M;
      IWORK = ITAUP + M;

      // Bidiagonalize A
      // (Workspace: need 3*M + N, prefer 3*M + (M + N)*NB)

      dgebrd(M, N, A, LDA, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP), WORK(IWORK),
          LWORK - IWORK + 1, IERR);
      if (WNTUAS) {
        // If left singular vectors desired in U, copy result to U
        // and generate left bidiagonalizing vectors in U
        // (Workspace: need 4*M-1, prefer 3*M + (M-1)*NB)

        dlacpy('L', M, M, A, LDA, U, LDU);
        dorgbr('Q', M, M, N, U, LDU, WORK(ITAUQ), WORK(IWORK),
            LWORK - IWORK + 1, IERR);
      }
      if (WNTVAS) {
        // If right singular vectors desired in VT, copy result to
        // VT and generate right bidiagonalizing vectors in VT
        // (Workspace: need 3*M + NRVT, prefer 3*M + NRVT*NB)

        dlacpy('U', M, N, A, LDA, VT, LDVT);
        if (WNTVA) NRVT = N;
        if (WNTVS) NRVT = M;
        dorgbr('P', NRVT, N, M, VT, LDVT, WORK(ITAUP), WORK(IWORK),
            LWORK - IWORK + 1, IERR);
      }
      if (WNTUO) {
        // If left singular vectors desired in A, generate left
        // bidiagonalizing vectors in A
        // (Workspace: need 4*M-1, prefer 3*M + (M-1)*NB)

        dorgbr('Q', M, M, N, A, LDA, WORK(ITAUQ), WORK(IWORK),
            LWORK - IWORK + 1, IERR);
      }
      if (WNTVO) {
        // If right singular vectors desired in A, generate right
        // bidiagonalizing vectors in A
        // (Workspace: need 4*M, prefer 3*M + M*NB)

        dorgbr('P', M, N, M, A, LDA, WORK(ITAUP), WORK(IWORK),
            LWORK - IWORK + 1, IERR);
      }
      IWORK = IE + M;
      if (WNTUAS || WNTUO) NRU = M;
      if (WNTUN) NRU = 0;
      if (WNTVAS || WNTVO) NCVT = N;
      if (WNTVN) NCVT = 0;
      if (!WNTUO && !WNTVO) {
        // Perform bidiagonal QR iteration, if desired, computing
        // left singular vectors in U and computing right singular
        // vectors in VT
        // (Workspace: need BDSPAC)

        dbdsqr('L', M, NCVT, NRU, 0, S, WORK(IE), VT, LDVT, U, LDU,
            DUM.asMatrix(1), 1, WORK(IWORK), INFO);
      } else if (!WNTUO && WNTVO) {
        // Perform bidiagonal QR iteration, if desired, computing
        // left singular vectors in U and computing right singular
        // vectors in A
        // (Workspace: need BDSPAC)

        dbdsqr('L', M, NCVT, NRU, 0, S, WORK(IE), A, LDA, U, LDU,
            DUM.asMatrix(1), 1, WORK(IWORK), INFO);
      } else {
        // Perform bidiagonal QR iteration, if desired, computing
        // left singular vectors in A and computing right singular
        // vectors in VT
        // (Workspace: need BDSPAC)

        dbdsqr('L', M, NCVT, NRU, 0, S, WORK(IE), VT, LDVT, A, LDA,
            DUM.asMatrix(1), 1, WORK(IWORK), INFO);
      }
    }
  }

  // If DBDSQR failed to converge, copy unconverged superdiagonals
  // to WORK( 2:MINMN )

  if (INFO.value != 0) {
    if (IE > 2) {
      for (I = 1; I <= MINMN - 1; I++) {
        WORK[I + 1] = WORK[I + IE - 1];
      }
    }
    if (IE < 2) {
      for (I = MINMN - 1; I >= 1; I--) {
        WORK[I + 1] = WORK[I + IE - 1];
      }
    }
  }

  // Undo scaling if necessary

  if (ISCL == 1) {
    if (ANRM > BIGNUM) {
      dlascl('G', 0, 0, BIGNUM, ANRM, MINMN, 1, S.asMatrix(MINMN), MINMN, IERR);
    }
    if (INFO.value != 0 && ANRM > BIGNUM) {
      dlascl('G', 0, 0, BIGNUM, ANRM, MINMN - 1, 1, WORK(2).asMatrix(MINMN),
          MINMN, IERR);
    }
    if (ANRM < SMLNUM) {
      dlascl('G', 0, 0, SMLNUM, ANRM, MINMN, 1, S.asMatrix(MINMN), MINMN, IERR);
    }
    if (INFO.value != 0 && ANRM < SMLNUM) {
      dlascl('G', 0, 0, SMLNUM, ANRM, MINMN - 1, 1, WORK(2).asMatrix(MINMN),
          MINMN, IERR);
    }
  }

  // Return optimal workspace in WORK(1)

  WORK[1] = MAXWRK.toDouble();
}
