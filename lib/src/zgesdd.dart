// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/dbdsdc.dart';
import 'package:dart_lapack/src/disnan.dart';
import 'package:dart_lapack/src/dlascl.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/droundup_lwork.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zgeqrf.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zgebrd.dart';
import 'package:dart_lapack/src/zgelqf.dart';
import 'package:dart_lapack/src/zlacp2.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlacrm.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlarcm.dart';
import 'package:dart_lapack/src/zlascl.dart';
import 'package:dart_lapack/src/zlaset.dart';
import 'package:dart_lapack/src/zungbr.dart';
import 'package:dart_lapack/src/zunglq.dart';
import 'package:dart_lapack/src/zungqr.dart';
import 'package:dart_lapack/src/zunmbr.dart';

void zgesdd(
  final String JOBZ,
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> S_,
  final Matrix<Complex> U_,
  final int LDU,
  final Matrix<Complex> VT_,
  final int LDVT,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final U = U_.having(ld: LDU);
  final VT = VT_.having(ld: LDVT);
  final S = S_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool LQUERY, WNTQA, WNTQAS, WNTQN, WNTQO, WNTQS;
  int BLK,
      CHUNK = 0,
      I,
      IE = 0,
      IL,
      IR,
      IRU,
      IRVT,
      ISCL,
      ITAU,
      ITAUP,
      ITAUQ,
      IU,
      IVT,
      LDWKVT,
      LDWRKL,
      LDWRKR,
      LDWRKU,
      MAXWRK,
      MINMN,
      MINWRK,
      MNTHR1,
      MNTHR2,
      NRWORK,
      NWORK,
      WRKBL;
  int LWORK_ZGEBRD_MN,
      LWORK_ZGEBRD_MM,
      LWORK_ZGEBRD_NN,
      LWORK_ZGELQF_MN,
      LWORK_ZGEQRF_MN,
      LWORK_ZUNGBR_P_MN,
      LWORK_ZUNGBR_P_NN,
      LWORK_ZUNGBR_Q_MN,
      LWORK_ZUNGBR_Q_MM,
      LWORK_ZUNGLQ_MN,
      LWORK_ZUNGLQ_NN,
      LWORK_ZUNGQR_MM,
      LWORK_ZUNGQR_MN,
      LWORK_ZUNMBR_PRC_MM,
      LWORK_ZUNMBR_QLN_MM,
      LWORK_ZUNMBR_PRC_MN,
      LWORK_ZUNMBR_QLN_MN,
      LWORK_ZUNMBR_PRC_NN,
      LWORK_ZUNMBR_QLN_NN;
  double ANRM, BIGNUM, EPS, SMLNUM;
  final IDUM = Array<int>(1);
  final DUM = Array<double>(1);
  final CDUM = Array<Complex>(1);
  final IERR = Box(0);

  // Test the input arguments

  INFO.value = 0;
  MINMN = min(M, N);
  MNTHR1 = MINMN * 17.0 ~/ 9.0;
  MNTHR2 = MINMN * 5.0 ~/ 3.0;
  WNTQA = lsame(JOBZ, 'A');
  WNTQS = lsame(JOBZ, 'S');
  WNTQAS = WNTQA || WNTQS;
  WNTQO = lsame(JOBZ, 'O');
  WNTQN = lsame(JOBZ, 'N');
  LQUERY = (LWORK == -1);
  MINWRK = 1;
  MAXWRK = 1;

  if (!(WNTQA || WNTQS || WNTQO || WNTQN)) {
    INFO.value = -1;
  } else if (M < 0) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, M)) {
    INFO.value = -5;
  } else if (LDU < 1 || (WNTQAS && LDU < M) || (WNTQO && M < N && LDU < M)) {
    INFO.value = -8;
  } else if (LDVT < 1 ||
      (WNTQA && LDVT < N) ||
      (WNTQS && LDVT < MINMN) ||
      (WNTQO && M >= N && LDVT < N)) {
    INFO.value = -10;
  }

  // Compute workspace
  //   Note: Comments in the code beginning "Workspace:" describe the
  //   minimal amount of workspace allocated at that point in the code,
  //   as well as the preferred amount for good performance.
  //   CWorkspace refers to complex workspace, and RWorkspace to
  //   real workspace. NB refers to the optimal block size for the
  //   immediately following subroutine, as returned by ILAENV.)

  if (INFO.value == 0) {
    MINWRK = 1;
    MAXWRK = 1;
    if (M >= N && MINMN > 0) {
      // There is no complex work space needed for bidiagonal SVD
      // The real work space needed for bidiagonal SVD (dbdsdc) is
      // BDSPAC = 3*N*N + 4*N for singular values and vectors;
      // BDSPAC = 4*N         for singular values only;
      // not including e, RU, and RVT matrices.

      // Compute space preferred for each routine
      zgebrd(M, N, CDUM(1).asMatrix(M), M, DUM(1), DUM(1), CDUM(1), CDUM(1),
          CDUM(1), -1, IERR);
      LWORK_ZGEBRD_MN = CDUM[1].toInt();

      zgebrd(N, N, CDUM(1).asMatrix(N), N, DUM(1), DUM(1), CDUM(1), CDUM(1),
          CDUM(1), -1, IERR);
      LWORK_ZGEBRD_NN = CDUM[1].toInt();

      zgeqrf(M, N, CDUM(1).asMatrix(M), M, CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZGEQRF_MN = CDUM[1].toInt();

      zungbr('P', N, N, N, CDUM(1).asMatrix(N), N, CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZUNGBR_P_NN = CDUM[1].toInt();

      zungbr('Q', M, M, N, CDUM(1).asMatrix(M), M, CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZUNGBR_Q_MM = CDUM[1].toInt();

      zungbr('Q', M, N, N, CDUM(1).asMatrix(M), M, CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZUNGBR_Q_MN = CDUM[1].toInt();

      zungqr(M, M, N, CDUM(1).asMatrix(M), M, CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZUNGQR_MM = CDUM[1].toInt();

      zungqr(M, N, N, CDUM(1).asMatrix(M), M, CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZUNGQR_MN = CDUM[1].toInt();

      zunmbr('P', 'R', 'C', N, N, N, CDUM(1).asMatrix(N), N, CDUM(1),
          CDUM(1).asMatrix(N), N, CDUM(1), -1, IERR);
      LWORK_ZUNMBR_PRC_NN = CDUM[1].toInt();

      zunmbr('Q', 'L', 'N', M, M, N, CDUM(1).asMatrix(M), M, CDUM(1),
          CDUM(1).asMatrix(M), M, CDUM(1), -1, IERR);
      LWORK_ZUNMBR_QLN_MM = CDUM[1].toInt();

      zunmbr('Q', 'L', 'N', M, N, N, CDUM(1).asMatrix(M), M, CDUM(1),
          CDUM(1).asMatrix(M), M, CDUM(1), -1, IERR);
      LWORK_ZUNMBR_QLN_MN = CDUM[1].toInt();

      zunmbr('Q', 'L', 'N', N, N, N, CDUM(1).asMatrix(N), N, CDUM(1),
          CDUM(1).asMatrix(N), N, CDUM(1), -1, IERR);
      LWORK_ZUNMBR_QLN_NN = CDUM[1].toInt();

      if (M >= MNTHR1) {
        if (WNTQN) {
          // Path 1 (M >> N, JOBZ='N')

          MAXWRK = N + LWORK_ZGEQRF_MN;
          MAXWRK = max(MAXWRK, 2 * N + LWORK_ZGEBRD_NN);
          MINWRK = 3 * N;
        } else if (WNTQO) {
          // Path 2 (M >> N, JOBZ='O')

          WRKBL = N + LWORK_ZGEQRF_MN;
          WRKBL = max(WRKBL, N + LWORK_ZUNGQR_MN);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZGEBRD_NN);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZUNMBR_QLN_NN);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZUNMBR_PRC_NN);
          MAXWRK = M * N + N * N + WRKBL;
          MINWRK = 2 * N * N + 3 * N;
        } else if (WNTQS) {
          // Path 3 (M >> N, JOBZ='S')

          WRKBL = N + LWORK_ZGEQRF_MN;
          WRKBL = max(WRKBL, N + LWORK_ZUNGQR_MN);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZGEBRD_NN);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZUNMBR_QLN_NN);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZUNMBR_PRC_NN);
          MAXWRK = N * N + WRKBL;
          MINWRK = N * N + 3 * N;
        } else if (WNTQA) {
          // Path 4 (M >> N, JOBZ='A')

          WRKBL = N + LWORK_ZGEQRF_MN;
          WRKBL = max(WRKBL, N + LWORK_ZUNGQR_MM);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZGEBRD_NN);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZUNMBR_QLN_NN);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZUNMBR_PRC_NN);
          MAXWRK = N * N + WRKBL;
          MINWRK = N * N + max(3 * N, N + M);
        }
      } else if (M >= MNTHR2) {
        // Path 5 (M >> N, but not as much as MNTHR1)

        MAXWRK = 2 * N + LWORK_ZGEBRD_MN;
        MINWRK = 2 * N + M;
        if (WNTQO) {
          // Path 5o (M >> N, JOBZ='O')
          MAXWRK = max(MAXWRK, 2 * N + LWORK_ZUNGBR_P_NN);
          MAXWRK = max(MAXWRK, 2 * N + LWORK_ZUNGBR_Q_MN);
          MAXWRK += M * N;
          MINWRK += N * N;
        } else if (WNTQS) {
          // Path 5s (M >> N, JOBZ='S')
          MAXWRK = max(MAXWRK, 2 * N + LWORK_ZUNGBR_P_NN);
          MAXWRK = max(MAXWRK, 2 * N + LWORK_ZUNGBR_Q_MN);
        } else if (WNTQA) {
          // Path 5a (M >> N, JOBZ='A')
          MAXWRK = max(MAXWRK, 2 * N + LWORK_ZUNGBR_P_NN);
          MAXWRK = max(MAXWRK, 2 * N + LWORK_ZUNGBR_Q_MM);
        }
      } else {
        // Path 6 (M >= N, but not much larger)

        MAXWRK = 2 * N + LWORK_ZGEBRD_MN;
        MINWRK = 2 * N + M;
        if (WNTQO) {
          // Path 6o (M >= N, JOBZ='O')
          MAXWRK = max(MAXWRK, 2 * N + LWORK_ZUNMBR_PRC_NN);
          MAXWRK = max(MAXWRK, 2 * N + LWORK_ZUNMBR_QLN_MN);
          MAXWRK += M * N;
          MINWRK += N * N;
        } else if (WNTQS) {
          // Path 6s (M >= N, JOBZ='S')
          MAXWRK = max(MAXWRK, 2 * N + LWORK_ZUNMBR_QLN_MN);
          MAXWRK = max(MAXWRK, 2 * N + LWORK_ZUNMBR_PRC_NN);
        } else if (WNTQA) {
          // Path 6a (M >= N, JOBZ='A')
          MAXWRK = max(MAXWRK, 2 * N + LWORK_ZUNMBR_QLN_MM);
          MAXWRK = max(MAXWRK, 2 * N + LWORK_ZUNMBR_PRC_NN);
        }
      }
    } else if (MINMN > 0) {
      // There is no complex work space needed for bidiagonal SVD
      // The real work space needed for bidiagonal SVD (dbdsdc) is
      // BDSPAC = 3*M*M + 4*M for singular values and vectors;
      // BDSPAC = 4*M         for singular values only;
      // not including e, RU, and RVT matrices.

      // Compute space preferred for each routine
      zgebrd(M, N, CDUM(1).asMatrix(M), M, DUM(1), DUM(1), CDUM(1), CDUM(1),
          CDUM(1), -1, IERR);
      LWORK_ZGEBRD_MN = CDUM[1].toInt();

      zgebrd(M, M, CDUM(1).asMatrix(M), M, DUM(1), DUM(1), CDUM(1), CDUM(1),
          CDUM(1), -1, IERR);
      LWORK_ZGEBRD_MM = CDUM[1].toInt();

      zgelqf(M, N, CDUM(1).asMatrix(M), M, CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZGELQF_MN = CDUM[1].toInt();

      zungbr('P', M, N, M, CDUM(1).asMatrix(M), M, CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZUNGBR_P_MN = CDUM[1].toInt();

      zungbr('P', N, N, M, CDUM(1).asMatrix(N), N, CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZUNGBR_P_NN = CDUM[1].toInt();

      zungbr('Q', M, M, N, CDUM(1).asMatrix(M), M, CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZUNGBR_Q_MM = CDUM[1].toInt();

      zunglq(M, N, M, CDUM(1).asMatrix(M), M, CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZUNGLQ_MN = CDUM[1].toInt();

      zunglq(N, N, M, CDUM(1).asMatrix(N), N, CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZUNGLQ_NN = CDUM[1].toInt();

      zunmbr('P', 'R', 'C', M, M, M, CDUM(1).asMatrix(M), M, CDUM(1),
          CDUM(1).asMatrix(M), M, CDUM(1), -1, IERR);
      LWORK_ZUNMBR_PRC_MM = CDUM[1].toInt();

      zunmbr('P', 'R', 'C', M, N, M, CDUM(1).asMatrix(M), M, CDUM(1),
          CDUM(1).asMatrix(M), M, CDUM(1), -1, IERR);
      LWORK_ZUNMBR_PRC_MN = CDUM[1].toInt();

      zunmbr('P', 'R', 'C', N, N, M, CDUM(1).asMatrix(N), N, CDUM(1),
          CDUM(1).asMatrix(N), N, CDUM(1), -1, IERR);
      LWORK_ZUNMBR_PRC_NN = CDUM[1].toInt();

      zunmbr('Q', 'L', 'N', M, M, M, CDUM(1).asMatrix(M), M, CDUM(1),
          CDUM(1).asMatrix(M), M, CDUM(1), -1, IERR);
      LWORK_ZUNMBR_QLN_MM = CDUM[1].toInt();

      if (N >= MNTHR1) {
        if (WNTQN) {
          // Path 1t (N >> M, JOBZ='N')

          MAXWRK = M + LWORK_ZGELQF_MN;
          MAXWRK = max(MAXWRK, 2 * M + LWORK_ZGEBRD_MM);
          MINWRK = 3 * M;
        } else if (WNTQO) {
          // Path 2t (N >> M, JOBZ='O')

          WRKBL = M + LWORK_ZGELQF_MN;
          WRKBL = max(WRKBL, M + LWORK_ZUNGLQ_MN);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZGEBRD_MM);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZUNMBR_QLN_MM);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZUNMBR_PRC_MM);
          MAXWRK = M * N + M * M + WRKBL;
          MINWRK = 2 * M * M + 3 * M;
        } else if (WNTQS) {
          // Path 3t (N >> M, JOBZ='S')

          WRKBL = M + LWORK_ZGELQF_MN;
          WRKBL = max(WRKBL, M + LWORK_ZUNGLQ_MN);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZGEBRD_MM);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZUNMBR_QLN_MM);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZUNMBR_PRC_MM);
          MAXWRK = M * M + WRKBL;
          MINWRK = M * M + 3 * M;
        } else if (WNTQA) {
          // Path 4t (N >> M, JOBZ='A')

          WRKBL = M + LWORK_ZGELQF_MN;
          WRKBL = max(WRKBL, M + LWORK_ZUNGLQ_NN);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZGEBRD_MM);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZUNMBR_QLN_MM);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZUNMBR_PRC_MM);
          MAXWRK = M * M + WRKBL;
          MINWRK = M * M + max(3 * M, M + N);
        }
      } else if (N >= MNTHR2) {
        // Path 5t (N >> M, but not as much as MNTHR1)

        MAXWRK = 2 * M + LWORK_ZGEBRD_MN;
        MINWRK = 2 * M + N;
        if (WNTQO) {
          // Path 5to (N >> M, JOBZ='O')
          MAXWRK = max(MAXWRK, 2 * M + LWORK_ZUNGBR_Q_MM);
          MAXWRK = max(MAXWRK, 2 * M + LWORK_ZUNGBR_P_MN);
          MAXWRK += M * N;
          MINWRK += M * M;
        } else if (WNTQS) {
          // Path 5ts (N >> M, JOBZ='S')
          MAXWRK = max(MAXWRK, 2 * M + LWORK_ZUNGBR_Q_MM);
          MAXWRK = max(MAXWRK, 2 * M + LWORK_ZUNGBR_P_MN);
        } else if (WNTQA) {
          // Path 5ta (N >> M, JOBZ='A')
          MAXWRK = max(MAXWRK, 2 * M + LWORK_ZUNGBR_Q_MM);
          MAXWRK = max(MAXWRK, 2 * M + LWORK_ZUNGBR_P_NN);
        }
      } else {
        // Path 6t (N > M, but not much larger)

        MAXWRK = 2 * M + LWORK_ZGEBRD_MN;
        MINWRK = 2 * M + N;
        if (WNTQO) {
          // Path 6to (N > M, JOBZ='O')
          MAXWRK = max(MAXWRK, 2 * M + LWORK_ZUNMBR_QLN_MM);
          MAXWRK = max(MAXWRK, 2 * M + LWORK_ZUNMBR_PRC_MN);
          MAXWRK += M * N;
          MINWRK += M * M;
        } else if (WNTQS) {
          // Path 6ts (N > M, JOBZ='S')
          MAXWRK = max(MAXWRK, 2 * M + LWORK_ZUNMBR_QLN_MM);
          MAXWRK = max(MAXWRK, 2 * M + LWORK_ZUNMBR_PRC_MN);
        } else if (WNTQA) {
          // Path 6ta (N > M, JOBZ='A')
          MAXWRK = max(MAXWRK, 2 * M + LWORK_ZUNMBR_QLN_MM);
          MAXWRK = max(MAXWRK, 2 * M + LWORK_ZUNMBR_PRC_NN);
        }
      }
    }
    MAXWRK = max(MAXWRK, MINWRK);
  }
  if (INFO.value == 0) {
    WORK[1] = droundup_lwork(MAXWRK).toComplex();
    if (LWORK < MINWRK && !LQUERY) {
      INFO.value = -12;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZGESDD', -INFO.value);
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

  ANRM = zlange('M', M, N, A, LDA, DUM);
  if (disnan(ANRM)) {
    INFO.value = -4;
    return;
  }
  ISCL = 0;
  if (ANRM > ZERO && ANRM < SMLNUM) {
    ISCL = 1;
    zlascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, IERR);
  } else if (ANRM > BIGNUM) {
    ISCL = 1;
    zlascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, IERR);
  }

  if (M >= N) {
    // A has at least as many rows as columns. If A has sufficiently
    // more rows than columns, first reduce using the QR
    // decomposition (if sufficient workspace available)

    if (M >= MNTHR1) {
      if (WNTQN) {
        // Path 1 (M >> N, JOBZ='N')
        // No singular vectors to be computed

        ITAU = 1;
        NWORK = ITAU + N;

        // Compute A=Q*R
        // CWorkspace: need   N [tau] + N    [work]
        // CWorkspace: prefer N [tau] + N*NB [work]
        // RWorkspace: need   0

        zgeqrf(M, N, A, LDA, WORK(ITAU), WORK(NWORK), LWORK - NWORK + 1, IERR);

        // Zero out below R

        zlaset('L', N - 1, N - 1, Complex.zero, Complex.zero, A(2, 1), LDA);
        IE = 1;
        ITAUQ = 1;
        ITAUP = ITAUQ + N;
        NWORK = ITAUP + N;

        // Bidiagonalize R in A
        // CWorkspace: need   2*N [tauq, taup] + N      [work]
        // CWorkspace: prefer 2*N [tauq, taup] + 2*N*NB [work]
        // RWorkspace: need   N [e]

        zgebrd(N, N, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
            WORK(NWORK), LWORK - NWORK + 1, IERR);
        NRWORK = IE + N;

        // Perform bidiagonal SVD, compute singular values only
        // CWorkspace: need   0
        // RWorkspace: need   N [e] + BDSPAC

        dbdsdc('U', 'N', N, S, RWORK(IE), DUM.asMatrix(1), 1, DUM.asMatrix(1),
            1, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);
      } else if (WNTQO) {
        // Path 2 (M >> N, JOBZ='O')
        // N left singular vectors to be overwritten on A and
        // N right singular vectors to be computed in VT

        IU = 1;

        // WORK(IU) is N by N

        LDWRKU = N;
        IR = IU + LDWRKU * N;
        if (LWORK >= M * N + N * N + 3 * N) {
          // WORK(IR) is M by N

          LDWRKR = M;
        } else {
          LDWRKR = (LWORK - N * N - 3 * N) ~/ N;
        }
        ITAU = IR + LDWRKR * N;
        NWORK = ITAU + N;

        // Compute A=Q*R
        // CWorkspace: need   N*N [U] + N*N [R] + N [tau] + N    [work]
        // CWorkspace: prefer N*N [U] + N*N [R] + N [tau] + N*NB [work]
        // RWorkspace: need   0

        zgeqrf(M, N, A, LDA, WORK(ITAU), WORK(NWORK), LWORK - NWORK + 1, IERR);

        // Copy R to WORK( IR ), zeroing out below it

        zlacpy('U', N, N, A, LDA, WORK(IR).asMatrix(LDWRKR), LDWRKR);
        zlaset('L', N - 1, N - 1, Complex.zero, Complex.zero,
            WORK(IR + 1).asMatrix(LDWRKR), LDWRKR);

        // Generate Q in A
        // CWorkspace: need   N*N [U] + N*N [R] + N [tau] + N    [work]
        // CWorkspace: prefer N*N [U] + N*N [R] + N [tau] + N*NB [work]
        // RWorkspace: need   0

        zungqr(
            M, N, N, A, LDA, WORK(ITAU), WORK(NWORK), LWORK - NWORK + 1, IERR);
        IE = 1;
        ITAUQ = ITAU;
        ITAUP = ITAUQ + N;
        NWORK = ITAUP + N;

        // Bidiagonalize R in WORK(IR)
        // CWorkspace: need   N*N [U] + N*N [R] + 2*N [tauq, taup] + N      [work]
        // CWorkspace: prefer N*N [U] + N*N [R] + 2*N [tauq, taup] + 2*N*NB [work]
        // RWorkspace: need   N [e]

        zgebrd(N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, S, RWORK(IE),
            WORK(ITAUQ), WORK(ITAUP), WORK(NWORK), LWORK - NWORK + 1, IERR);

        // Perform bidiagonal SVD, computing left singular vectors
        // of R in WORK(IRU) and computing right singular vectors
        // of R in WORK(IRVT)
        // CWorkspace: need   0
        // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC

        IRU = IE + N;
        IRVT = IRU + N * N;
        NRWORK = IRVT + N * N;
        dbdsdc('U', 'I', N, S, RWORK(IE), RWORK(IRU).asMatrix(N), N,
            RWORK(IRVT).asMatrix(N), N, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);

        // Copy real matrix RWORK(IRU) to complex matrix WORK(IU)
        // Overwrite WORK(IU) by the left singular vectors of R
        // CWorkspace: need   N*N [U] + N*N [R] + 2*N [tauq, taup] + N    [work]
        // CWorkspace: prefer N*N [U] + N*N [R] + 2*N [tauq, taup] + N*NB [work]
        // RWorkspace: need   0

        zlacp2('F', N, N, RWORK(IRU).asMatrix(N), N, WORK(IU).asMatrix(LDWRKU),
            LDWRKU);
        zunmbr(
            'Q',
            'L',
            'N',
            N,
            N,
            N,
            WORK(IR).asMatrix(LDWRKR),
            LDWRKR,
            WORK(ITAUQ),
            WORK(IU).asMatrix(LDWRKU),
            LDWRKU,
            WORK(NWORK),
            LWORK - NWORK + 1,
            IERR);

        // Copy real matrix RWORK(IRVT) to complex matrix VT
        // Overwrite VT by the right singular vectors of R
        // CWorkspace: need   N*N [U] + N*N [R] + 2*N [tauq, taup] + N    [work]
        // CWorkspace: prefer N*N [U] + N*N [R] + 2*N [tauq, taup] + N*NB [work]
        // RWorkspace: need   0

        zlacp2('F', N, N, RWORK(IRVT).asMatrix(N), N, VT, LDVT);
        zunmbr('P', 'R', 'C', N, N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR,
            WORK(ITAUP), VT, LDVT, WORK(NWORK), LWORK - NWORK + 1, IERR);

        // Multiply Q in A by left singular vectors of R in
        // WORK(IU), storing result in WORK(IR) and copying to A
        // CWorkspace: need   N*N [U] + N*N [R]
        // CWorkspace: prefer N*N [U] + M*N [R]
        // RWorkspace: need   0

        for (I = 1; I <= M; I += LDWRKR) {
          CHUNK = min(M - I + 1, LDWRKR);
          zgemm(
              'N',
              'N',
              CHUNK,
              N,
              N,
              Complex.one,
              A(I, 1),
              LDA,
              WORK(IU).asMatrix(LDWRKU),
              LDWRKU,
              Complex.zero,
              WORK(IR).asMatrix(LDWRKR),
              LDWRKR);
          zlacpy(
              'F', CHUNK, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, A(I, 1), LDA);
        }
      } else if (WNTQS) {
        // Path 3 (M >> N, JOBZ='S')
        // N left singular vectors to be computed in U and
        // N right singular vectors to be computed in VT

        IR = 1;

        // WORK(IR) is N by N

        LDWRKR = N;
        ITAU = IR + LDWRKR * N;
        NWORK = ITAU + N;

        // Compute A=Q*R
        // CWorkspace: need   N*N [R] + N [tau] + N    [work]
        // CWorkspace: prefer N*N [R] + N [tau] + N*NB [work]
        // RWorkspace: need   0

        zgeqrf(M, N, A, LDA, WORK(ITAU), WORK(NWORK), LWORK - NWORK + 1, IERR);

        // Copy R to WORK(IR), zeroing out below it

        zlacpy('U', N, N, A, LDA, WORK(IR).asMatrix(LDWRKR), LDWRKR);
        zlaset('L', N - 1, N - 1, Complex.zero, Complex.zero,
            WORK(IR + 1).asMatrix(LDWRKR), LDWRKR);

        // Generate Q in A
        // CWorkspace: need   N*N [R] + N [tau] + N    [work]
        // CWorkspace: prefer N*N [R] + N [tau] + N*NB [work]
        // RWorkspace: need   0

        zungqr(
            M, N, N, A, LDA, WORK(ITAU), WORK(NWORK), LWORK - NWORK + 1, IERR);
        IE = 1;
        ITAUQ = ITAU;
        ITAUP = ITAUQ + N;
        NWORK = ITAUP + N;

        // Bidiagonalize R in WORK(IR)
        // CWorkspace: need   N*N [R] + 2*N [tauq, taup] + N      [work]
        // CWorkspace: prefer N*N [R] + 2*N [tauq, taup] + 2*N*NB [work]
        // RWorkspace: need   N [e]

        zgebrd(N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, S, RWORK(IE),
            WORK(ITAUQ), WORK(ITAUP), WORK(NWORK), LWORK - NWORK + 1, IERR);

        // Perform bidiagonal SVD, computing left singular vectors
        // of bidiagonal matrix in RWORK(IRU) and computing right
        // singular vectors of bidiagonal matrix in RWORK(IRVT)
        // CWorkspace: need   0
        // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC

        IRU = IE + N;
        IRVT = IRU + N * N;
        NRWORK = IRVT + N * N;
        dbdsdc('U', 'I', N, S, RWORK(IE), RWORK(IRU).asMatrix(N), N,
            RWORK(IRVT).asMatrix(N), N, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);

        // Copy real matrix RWORK(IRU) to complex matrix U
        // Overwrite U by left singular vectors of R
        // CWorkspace: need   N*N [R] + 2*N [tauq, taup] + N    [work]
        // CWorkspace: prefer N*N [R] + 2*N [tauq, taup] + N*NB [work]
        // RWorkspace: need   0

        zlacp2('F', N, N, RWORK(IRU).asMatrix(N), N, U, LDU);
        zunmbr('Q', 'L', 'N', N, N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR,
            WORK(ITAUQ), U, LDU, WORK(NWORK), LWORK - NWORK + 1, IERR);

        // Copy real matrix RWORK(IRVT) to complex matrix VT
        // Overwrite VT by right singular vectors of R
        // CWorkspace: need   N*N [R] + 2*N [tauq, taup] + N    [work]
        // CWorkspace: prefer N*N [R] + 2*N [tauq, taup] + N*NB [work]
        // RWorkspace: need   0

        zlacp2('F', N, N, RWORK(IRVT).asMatrix(N), N, VT, LDVT);
        zunmbr('P', 'R', 'C', N, N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR,
            WORK(ITAUP), VT, LDVT, WORK(NWORK), LWORK - NWORK + 1, IERR);

        // Multiply Q in A by left singular vectors of R in
        // WORK(IR), storing result in U
        // CWorkspace: need   N*N [R]
        // RWorkspace: need   0

        zlacpy('F', N, N, U, LDU, WORK(IR).asMatrix(LDWRKR), LDWRKR);
        zgemm('N', 'N', M, N, N, Complex.one, A, LDA, WORK(IR).asMatrix(LDWRKR),
            LDWRKR, Complex.zero, U, LDU);
      } else if (WNTQA) {
        // Path 4 (M >> N, JOBZ='A')
        // M left singular vectors to be computed in U and
        // N right singular vectors to be computed in VT

        IU = 1;

        // WORK(IU) is N by N

        LDWRKU = N;
        ITAU = IU + LDWRKU * N;
        NWORK = ITAU + N;

        // Compute A=Q*R, copying result to U
        // CWorkspace: need   N*N [U] + N [tau] + N    [work]
        // CWorkspace: prefer N*N [U] + N [tau] + N*NB [work]
        // RWorkspace: need   0

        zgeqrf(M, N, A, LDA, WORK(ITAU), WORK(NWORK), LWORK - NWORK + 1, IERR);
        zlacpy('L', M, N, A, LDA, U, LDU);

        // Generate Q in U
        // CWorkspace: need   N*N [U] + N [tau] + M    [work]
        // CWorkspace: prefer N*N [U] + N [tau] + M*NB [work]
        // RWorkspace: need   0

        zungqr(
            M, M, N, U, LDU, WORK(ITAU), WORK(NWORK), LWORK - NWORK + 1, IERR);

        // Produce R in A, zeroing out below it

        zlaset('L', N - 1, N - 1, Complex.zero, Complex.zero, A(2, 1), LDA);
        IE = 1;
        ITAUQ = ITAU;
        ITAUP = ITAUQ + N;
        NWORK = ITAUP + N;

        // Bidiagonalize R in A
        // CWorkspace: need   N*N [U] + 2*N [tauq, taup] + N      [work]
        // CWorkspace: prefer N*N [U] + 2*N [tauq, taup] + 2*N*NB [work]
        // RWorkspace: need   N [e]

        zgebrd(N, N, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
            WORK(NWORK), LWORK - NWORK + 1, IERR);
        IRU = IE + N;
        IRVT = IRU + N * N;
        NRWORK = IRVT + N * N;

        // Perform bidiagonal SVD, computing left singular vectors
        // of bidiagonal matrix in RWORK(IRU) and computing right
        // singular vectors of bidiagonal matrix in RWORK(IRVT)
        // CWorkspace: need   0
        // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC

        dbdsdc('U', 'I', N, S, RWORK(IE), RWORK(IRU).asMatrix(N), N,
            RWORK(IRVT).asMatrix(N), N, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);

        // Copy real matrix RWORK(IRU) to complex matrix WORK(IU)
        // Overwrite WORK(IU) by left singular vectors of R
        // CWorkspace: need   N*N [U] + 2*N [tauq, taup] + N    [work]
        // CWorkspace: prefer N*N [U] + 2*N [tauq, taup] + N*NB [work]
        // RWorkspace: need   0

        zlacp2('F', N, N, RWORK(IRU).asMatrix(N), N, WORK(IU).asMatrix(LDWRKU),
            LDWRKU);
        zunmbr(
            'Q',
            'L',
            'N',
            N,
            N,
            N,
            A,
            LDA,
            WORK(ITAUQ),
            WORK(IU).asMatrix(LDWRKU),
            LDWRKU,
            WORK(NWORK),
            LWORK - NWORK + 1,
            IERR);

        // Copy real matrix RWORK(IRVT) to complex matrix VT
        // Overwrite VT by right singular vectors of R
        // CWorkspace: need   N*N [U] + 2*N [tauq, taup] + N    [work]
        // CWorkspace: prefer N*N [U] + 2*N [tauq, taup] + N*NB [work]
        // RWorkspace: need   0

        zlacp2('F', N, N, RWORK(IRVT).asMatrix(N), N, VT, LDVT);
        zunmbr('P', 'R', 'C', N, N, N, A, LDA, WORK(ITAUP), VT, LDVT,
            WORK(NWORK), LWORK - NWORK + 1, IERR);

        // Multiply Q in U by left singular vectors of R in
        // WORK(IU), storing result in A
        // CWorkspace: need   N*N [U]
        // RWorkspace: need   0

        zgemm('N', 'N', M, N, N, Complex.one, U, LDU, WORK(IU).asMatrix(LDWRKU),
            LDWRKU, Complex.zero, A, LDA);

        // Copy left singular vectors of A from A to U

        zlacpy('F', M, N, A, LDA, U, LDU);
      }
    } else if (M >= MNTHR2) {
      // MNTHR2 <= M < MNTHR1

      // Path 5 (M >> N, but not as much as MNTHR1)
      // Reduce to bidiagonal form without QR decomposition, use
      // ZUNGBR and matrix multiplication to compute singular vectors

      IE = 1;
      NRWORK = IE + N;
      ITAUQ = 1;
      ITAUP = ITAUQ + N;
      NWORK = ITAUP + N;

      // Bidiagonalize A
      // CWorkspace: need   2*N [tauq, taup] + M        [work]
      // CWorkspace: prefer 2*N [tauq, taup] + (M+N)*NB [work]
      // RWorkspace: need   N [e]

      zgebrd(M, N, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP), WORK(NWORK),
          LWORK - NWORK + 1, IERR);
      if (WNTQN) {
        // Path 5n (M >> N, JOBZ='N')
        // Compute singular values only
        // CWorkspace: need   0
        // RWorkspace: need   N [e] + BDSPAC

        dbdsdc('U', 'N', N, S, RWORK(IE), DUM.asMatrix(1), 1, DUM.asMatrix(1),
            1, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);
      } else if (WNTQO) {
        IU = NWORK;
        IRU = NRWORK;
        IRVT = IRU + N * N;
        NRWORK = IRVT + N * N;

        // Path 5o (M >> N, JOBZ='O')
        // Copy A to VT, generate P**H
        // CWorkspace: need   2*N [tauq, taup] + N    [work]
        // CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
        // RWorkspace: need   0

        zlacpy('U', N, N, A, LDA, VT, LDVT);
        zungbr('P', N, N, N, VT, LDVT, WORK(ITAUP), WORK(NWORK),
            LWORK - NWORK + 1, IERR);

        // Generate Q in A
        // CWorkspace: need   2*N [tauq, taup] + N    [work]
        // CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
        // RWorkspace: need   0

        zungbr('Q', M, N, N, A, LDA, WORK(ITAUQ), WORK(NWORK),
            LWORK - NWORK + 1, IERR);

        if (LWORK >= M * N + 3 * N) {
          // WORK( IU ) is M by N

          LDWRKU = M;
        } else {
          // WORK(IU) is LDWRKU by N

          LDWRKU = (LWORK - 3 * N) ~/ N;
        }
        NWORK = IU + LDWRKU * N;

        // Perform bidiagonal SVD, computing left singular vectors
        // of bidiagonal matrix in RWORK(IRU) and computing right
        // singular vectors of bidiagonal matrix in RWORK(IRVT)
        // CWorkspace: need   0
        // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC

        dbdsdc('U', 'I', N, S, RWORK(IE), RWORK(IRU).asMatrix(N), N,
            RWORK(IRVT).asMatrix(N), N, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);

        // Multiply real matrix RWORK(IRVT) by P**H in VT,
        // storing the result in WORK(IU), copying to VT
        // CWorkspace: need   2*N [tauq, taup] + N*N [U]
        // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + 2*N*N [rwork]

        zlarcm(N, N, RWORK(IRVT).asMatrix(N), N, VT, LDVT,
            WORK(IU).asMatrix(LDWRKU), LDWRKU, RWORK(NRWORK));
        zlacpy('F', N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, VT, LDVT);

        // Multiply Q in A by real matrix RWORK(IRU), storing the
        // result in WORK(IU), copying to A
        // CWorkspace: need   2*N [tauq, taup] + N*N [U]
        // CWorkspace: prefer 2*N [tauq, taup] + M*N [U]
        // RWorkspace: need   N [e] + N*N [RU] + 2*N*N [rwork]
        // RWorkspace: prefer N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here

        NRWORK = IRVT;
        for (I = 1; I <= M; I += LDWRKU) {
          CHUNK = min(M - I + 1, LDWRKU);
          zlacrm(CHUNK, N, A(I, 1), LDA, RWORK(IRU).asMatrix(N), N,
              WORK(IU).asMatrix(LDWRKU), LDWRKU, RWORK(NRWORK));
          zlacpy(
              'F', CHUNK, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, A(I, 1), LDA);
        }
      } else if (WNTQS) {
        // Path 5s (M >> N, JOBZ='S')
        // Copy A to VT, generate P**H
        // CWorkspace: need   2*N [tauq, taup] + N    [work]
        // CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
        // RWorkspace: need   0

        zlacpy('U', N, N, A, LDA, VT, LDVT);
        zungbr('P', N, N, N, VT, LDVT, WORK(ITAUP), WORK(NWORK),
            LWORK - NWORK + 1, IERR);

        // Copy A to U, generate Q
        // CWorkspace: need   2*N [tauq, taup] + N    [work]
        // CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
        // RWorkspace: need   0

        zlacpy('L', M, N, A, LDA, U, LDU);
        zungbr('Q', M, N, N, U, LDU, WORK(ITAUQ), WORK(NWORK),
            LWORK - NWORK + 1, IERR);

        // Perform bidiagonal SVD, computing left singular vectors
        // of bidiagonal matrix in RWORK(IRU) and computing right
        // singular vectors of bidiagonal matrix in RWORK(IRVT)
        // CWorkspace: need   0
        // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC

        IRU = NRWORK;
        IRVT = IRU + N * N;
        NRWORK = IRVT + N * N;
        dbdsdc('U', 'I', N, S, RWORK(IE), RWORK(IRU).asMatrix(N), N,
            RWORK(IRVT).asMatrix(N), N, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);

        // Multiply real matrix RWORK(IRVT) by P**H in VT,
        // storing the result in A, copying to VT
        // CWorkspace: need   0
        // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + 2*N*N [rwork]

        zlarcm(
            N, N, RWORK(IRVT).asMatrix(N), N, VT, LDVT, A, LDA, RWORK(NRWORK));
        zlacpy('F', N, N, A, LDA, VT, LDVT);

        // Multiply Q in U by real matrix RWORK(IRU), storing the
        // result in A, copying to U
        // CWorkspace: need   0
        // RWorkspace: need   N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here

        NRWORK = IRVT;
        zlacrm(M, N, U, LDU, RWORK(IRU).asMatrix(N), N, A, LDA, RWORK(NRWORK));
        zlacpy('F', M, N, A, LDA, U, LDU);
      } else {
        // Path 5a (M >> N, JOBZ='A')
        // Copy A to VT, generate P**H
        // CWorkspace: need   2*N [tauq, taup] + N    [work]
        // CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
        // RWorkspace: need   0

        zlacpy('U', N, N, A, LDA, VT, LDVT);
        zungbr('P', N, N, N, VT, LDVT, WORK(ITAUP), WORK(NWORK),
            LWORK - NWORK + 1, IERR);

        // Copy A to U, generate Q
        // CWorkspace: need   2*N [tauq, taup] + M    [work]
        // CWorkspace: prefer 2*N [tauq, taup] + M*NB [work]
        // RWorkspace: need   0

        zlacpy('L', M, N, A, LDA, U, LDU);
        zungbr('Q', M, M, N, U, LDU, WORK(ITAUQ), WORK(NWORK),
            LWORK - NWORK + 1, IERR);

        // Perform bidiagonal SVD, computing left singular vectors
        // of bidiagonal matrix in RWORK(IRU) and computing right
        // singular vectors of bidiagonal matrix in RWORK(IRVT)
        // CWorkspace: need   0
        // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC

        IRU = NRWORK;
        IRVT = IRU + N * N;
        NRWORK = IRVT + N * N;
        dbdsdc('U', 'I', N, S, RWORK(IE), RWORK(IRU).asMatrix(N), N,
            RWORK(IRVT).asMatrix(N), N, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);

        // Multiply real matrix RWORK(IRVT) by P**H in VT,
        // storing the result in A, copying to VT
        // CWorkspace: need   0
        // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + 2*N*N [rwork]

        zlarcm(
            N, N, RWORK(IRVT).asMatrix(N), N, VT, LDVT, A, LDA, RWORK(NRWORK));
        zlacpy('F', N, N, A, LDA, VT, LDVT);

        // Multiply Q in U by real matrix RWORK(IRU), storing the
        // result in A, copying to U
        // CWorkspace: need   0
        // RWorkspace: need   N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here

        NRWORK = IRVT;
        zlacrm(M, N, U, LDU, RWORK(IRU).asMatrix(N), N, A, LDA, RWORK(NRWORK));
        zlacpy('F', M, N, A, LDA, U, LDU);
      }
    } else {
      // M < MNTHR2

      // Path 6 (M >= N, but not much larger)
      // Reduce to bidiagonal form without QR decomposition
      // Use ZUNMBR to compute singular vectors

      IE = 1;
      NRWORK = IE + N;
      ITAUQ = 1;
      ITAUP = ITAUQ + N;
      NWORK = ITAUP + N;

      // Bidiagonalize A
      // CWorkspace: need   2*N [tauq, taup] + M        [work]
      // CWorkspace: prefer 2*N [tauq, taup] + (M+N)*NB [work]
      // RWorkspace: need   N [e]

      zgebrd(M, N, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP), WORK(NWORK),
          LWORK - NWORK + 1, IERR);
      if (WNTQN) {
        // Path 6n (M >= N, JOBZ='N')
        // Compute singular values only
        // CWorkspace: need   0
        // RWorkspace: need   N [e] + BDSPAC

        dbdsdc('U', 'N', N, S, RWORK(IE), DUM.asMatrix(1), 1, DUM.asMatrix(1),
            1, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);
      } else if (WNTQO) {
        IU = NWORK;
        IRU = NRWORK;
        IRVT = IRU + N * N;
        NRWORK = IRVT + N * N;
        if (LWORK >= M * N + 3 * N) {
          // WORK( IU ) is M by N

          LDWRKU = M;
        } else {
          // WORK( IU ) is LDWRKU by N

          LDWRKU = (LWORK - 3 * N) ~/ N;
        }
        NWORK = IU + LDWRKU * N;

        // Path 6o (M >= N, JOBZ='O')
        // Perform bidiagonal SVD, computing left singular vectors
        // of bidiagonal matrix in RWORK(IRU) and computing right
        // singular vectors of bidiagonal matrix in RWORK(IRVT)
        // CWorkspace: need   0
        // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC

        dbdsdc('U', 'I', N, S, RWORK(IE), RWORK(IRU).asMatrix(N), N,
            RWORK(IRVT).asMatrix(N), N, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);

        // Copy real matrix RWORK(IRVT) to complex matrix VT
        // Overwrite VT by right singular vectors of A
        // CWorkspace: need   2*N [tauq, taup] + N*N [U] + N    [work]
        // CWorkspace: prefer 2*N [tauq, taup] + N*N [U] + N*NB [work]
        // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT]

        zlacp2('F', N, N, RWORK(IRVT).asMatrix(N), N, VT, LDVT);
        zunmbr('P', 'R', 'C', N, N, N, A, LDA, WORK(ITAUP), VT, LDVT,
            WORK(NWORK), LWORK - NWORK + 1, IERR);

        if (LWORK >= M * N + 3 * N) {
          // Path 6o-fast
          // Copy real matrix RWORK(IRU) to complex matrix WORK(IU)
          // Overwrite WORK(IU) by left singular vectors of A, copying
          // to A
          // CWorkspace: need   2*N [tauq, taup] + M*N [U] + N    [work]
          // CWorkspace: prefer 2*N [tauq, taup] + M*N [U] + N*NB [work]
          // RWorkspace: need   N [e] + N*N [RU]

          zlaset('F', M, N, Complex.zero, Complex.zero,
              WORK(IU).asMatrix(LDWRKU), LDWRKU);
          zlacp2('F', N, N, RWORK(IRU).asMatrix(N), N,
              WORK(IU).asMatrix(LDWRKU), LDWRKU);
          zunmbr(
              'Q',
              'L',
              'N',
              M,
              N,
              N,
              A,
              LDA,
              WORK(ITAUQ),
              WORK(IU).asMatrix(LDWRKU),
              LDWRKU,
              WORK(NWORK),
              LWORK - NWORK + 1,
              IERR);
          zlacpy('F', M, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, A, LDA);
        } else {
          // Path 6o-slow
          // Generate Q in A
          // CWorkspace: need   2*N [tauq, taup] + N*N [U] + N    [work]
          // CWorkspace: prefer 2*N [tauq, taup] + N*N [U] + N*NB [work]
          // RWorkspace: need   0

          zungbr('Q', M, N, N, A, LDA, WORK(ITAUQ), WORK(NWORK),
              LWORK - NWORK + 1, IERR);

          // Multiply Q in A by real matrix RWORK(IRU), storing the
          // result in WORK(IU), copying to A
          // CWorkspace: need   2*N [tauq, taup] + N*N [U]
          // CWorkspace: prefer 2*N [tauq, taup] + M*N [U]
          // RWorkspace: need   N [e] + N*N [RU] + 2*N*N [rwork]
          // RWorkspace: prefer N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here

          NRWORK = IRVT;
          for (I = 1; I <= M; I += LDWRKU) {
            CHUNK = min(M - I + 1, LDWRKU);
            zlacrm(CHUNK, N, A(I, 1), LDA, RWORK(IRU).asMatrix(N), N,
                WORK(IU).asMatrix(LDWRKU), LDWRKU, RWORK(NRWORK));
            zlacpy(
                'F', CHUNK, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, A(I, 1), LDA);
          }
        }
      } else if (WNTQS) {
        // Path 6s (M >= N, JOBZ='S')
        // Perform bidiagonal SVD, computing left singular vectors
        // of bidiagonal matrix in RWORK(IRU) and computing right
        // singular vectors of bidiagonal matrix in RWORK(IRVT)
        // CWorkspace: need   0
        // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC

        IRU = NRWORK;
        IRVT = IRU + N * N;
        NRWORK = IRVT + N * N;
        dbdsdc('U', 'I', N, S, RWORK(IE), RWORK(IRU).asMatrix(N), N,
            RWORK(IRVT).asMatrix(N), N, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);

        // Copy real matrix RWORK(IRU) to complex matrix U
        // Overwrite U by left singular vectors of A
        // CWorkspace: need   2*N [tauq, taup] + N    [work]
        // CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
        // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT]

        zlaset('F', M, N, Complex.zero, Complex.zero, U, LDU);
        zlacp2('F', N, N, RWORK(IRU).asMatrix(N), N, U, LDU);
        zunmbr('Q', 'L', 'N', M, N, N, A, LDA, WORK(ITAUQ), U, LDU, WORK(NWORK),
            LWORK - NWORK + 1, IERR);

        // Copy real matrix RWORK(IRVT) to complex matrix VT
        // Overwrite VT by right singular vectors of A
        // CWorkspace: need   2*N [tauq, taup] + N    [work]
        // CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
        // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT]

        zlacp2('F', N, N, RWORK(IRVT).asMatrix(N), N, VT, LDVT);
        zunmbr('P', 'R', 'C', N, N, N, A, LDA, WORK(ITAUP), VT, LDVT,
            WORK(NWORK), LWORK - NWORK + 1, IERR);
      } else {
        // Path 6a (M >= N, JOBZ='A')
        // Perform bidiagonal SVD, computing left singular vectors
        // of bidiagonal matrix in RWORK(IRU) and computing right
        // singular vectors of bidiagonal matrix in RWORK(IRVT)
        // CWorkspace: need   0
        // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC

        IRU = NRWORK;
        IRVT = IRU + N * N;
        NRWORK = IRVT + N * N;
        dbdsdc('U', 'I', N, S, RWORK(IE), RWORK(IRU).asMatrix(N), N,
            RWORK(IRVT).asMatrix(N), N, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);

        // Set the right corner of U to identity matrix

        zlaset('F', M, M, Complex.zero, Complex.zero, U, LDU);
        if (M > N) {
          zlaset('F', M - N, M - N, Complex.zero, Complex.one, U(N + 1, N + 1),
              LDU);
        }

        // Copy real matrix RWORK(IRU) to complex matrix U
        // Overwrite U by left singular vectors of A
        // CWorkspace: need   2*N [tauq, taup] + M    [work]
        // CWorkspace: prefer 2*N [tauq, taup] + M*NB [work]
        // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT]

        zlacp2('F', N, N, RWORK(IRU).asMatrix(N), N, U, LDU);
        zunmbr('Q', 'L', 'N', M, M, N, A, LDA, WORK(ITAUQ), U, LDU, WORK(NWORK),
            LWORK - NWORK + 1, IERR);

        // Copy real matrix RWORK(IRVT) to complex matrix VT
        // Overwrite VT by right singular vectors of A
        // CWorkspace: need   2*N [tauq, taup] + N    [work]
        // CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
        // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT]

        zlacp2('F', N, N, RWORK(IRVT).asMatrix(N), N, VT, LDVT);
        zunmbr('P', 'R', 'C', N, N, N, A, LDA, WORK(ITAUP), VT, LDVT,
            WORK(NWORK), LWORK - NWORK + 1, IERR);
      }
    }
  } else {
    // A has more columns than rows. If A has sufficiently more
    // columns than rows, first reduce using the LQ decomposition (if
    // sufficient workspace available)

    if (N >= MNTHR1) {
      if (WNTQN) {
        // Path 1t (N >> M, JOBZ='N')
        // No singular vectors to be computed

        ITAU = 1;
        NWORK = ITAU + M;

        // Compute A=L*Q
        // CWorkspace: need   M [tau] + M    [work]
        // CWorkspace: prefer M [tau] + M*NB [work]
        // RWorkspace: need   0

        zgelqf(M, N, A, LDA, WORK(ITAU), WORK(NWORK), LWORK - NWORK + 1, IERR);

        // Zero out above L

        zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero, A(1, 2), LDA);
        IE = 1;
        ITAUQ = 1;
        ITAUP = ITAUQ + M;
        NWORK = ITAUP + M;

        // Bidiagonalize L in A
        // CWorkspace: need   2*M [tauq, taup] + M      [work]
        // CWorkspace: prefer 2*M [tauq, taup] + 2*M*NB [work]
        // RWorkspace: need   M [e]

        zgebrd(M, M, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
            WORK(NWORK), LWORK - NWORK + 1, IERR);
        NRWORK = IE + M;

        // Perform bidiagonal SVD, compute singular values only
        // CWorkspace: need   0
        // RWorkspace: need   M [e] + BDSPAC

        dbdsdc('U', 'N', M, S, RWORK(IE), DUM.asMatrix(1), 1, DUM.asMatrix(1),
            1, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);
      } else if (WNTQO) {
        // Path 2t (N >> M, JOBZ='O')
        // M right singular vectors to be overwritten on A and
        // M left singular vectors to be computed in U

        IVT = 1;
        LDWKVT = M;

        // WORK(IVT) is M by M

        IL = IVT + LDWKVT * M;
        if (LWORK >= M * N + M * M + 3 * M) {
          // WORK(IL) M by N

          LDWRKL = M;
          CHUNK = N;
        } else {
          // WORK(IL) is M by CHUNK

          LDWRKL = M;
          CHUNK = (LWORK - M * M - 3 * M) ~/ M;
        }
        ITAU = IL + LDWRKL * CHUNK;
        NWORK = ITAU + M;

        // Compute A=L*Q
        // CWorkspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work]
        // CWorkspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work]
        // RWorkspace: need   0

        zgelqf(M, N, A, LDA, WORK(ITAU), WORK(NWORK), LWORK - NWORK + 1, IERR);

        // Copy L to WORK(IL), zeroing about above it

        zlacpy('L', M, M, A, LDA, WORK(IL).asMatrix(LDWRKL), LDWRKL);
        zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero,
            WORK(IL + LDWRKL).asMatrix(LDWRKL), LDWRKL);

        // Generate Q in A
        // CWorkspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work]
        // CWorkspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work]
        // RWorkspace: need   0

        zunglq(
            M, N, M, A, LDA, WORK(ITAU), WORK(NWORK), LWORK - NWORK + 1, IERR);
        IE = 1;
        ITAUQ = ITAU;
        ITAUP = ITAUQ + M;
        NWORK = ITAUP + M;

        // Bidiagonalize L in WORK(IL)
        // CWorkspace: need   M*M [VT] + M*M [L] + 2*M [tauq, taup] + M      [work]
        // CWorkspace: prefer M*M [VT] + M*M [L] + 2*M [tauq, taup] + 2*M*NB [work]
        // RWorkspace: need   M [e]

        zgebrd(M, M, WORK(IL).asMatrix(LDWRKL), LDWRKL, S, RWORK(IE),
            WORK(ITAUQ), WORK(ITAUP), WORK(NWORK), LWORK - NWORK + 1, IERR);

        // Perform bidiagonal SVD, computing left singular vectors
        // of bidiagonal matrix in RWORK(IRU) and computing right
        // singular vectors of bidiagonal matrix in RWORK(IRVT)
        // CWorkspace: need   0
        // RWorkspace: need   M [e] + M*M [RU] + M*M [RVT] + BDSPAC

        IRU = IE + M;
        IRVT = IRU + M * M;
        NRWORK = IRVT + M * M;
        dbdsdc('U', 'I', M, S, RWORK(IE), RWORK(IRU).asMatrix(M), M,
            RWORK(IRVT).asMatrix(M), M, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);

        // Copy real matrix RWORK(IRU) to complex matrix WORK(IU)
        // Overwrite WORK(IU) by the left singular vectors of L
        // CWorkspace: need   M*M [VT] + M*M [L] + 2*M [tauq, taup] + M    [work]
        // CWorkspace: prefer M*M [VT] + M*M [L] + 2*M [tauq, taup] + M*NB [work]
        // RWorkspace: need   0

        zlacp2('F', M, M, RWORK(IRU).asMatrix(M), M, U, LDU);
        zunmbr('Q', 'L', 'N', M, M, M, WORK(IL).asMatrix(LDWRKL), LDWRKL,
            WORK(ITAUQ), U, LDU, WORK(NWORK), LWORK - NWORK + 1, IERR);

        // Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT)
        // Overwrite WORK(IVT) by the right singular vectors of L
        // CWorkspace: need   M*M [VT] + M*M [L] + 2*M [tauq, taup] + M    [work]
        // CWorkspace: prefer M*M [VT] + M*M [L] + 2*M [tauq, taup] + M*NB [work]
        // RWorkspace: need   0

        zlacp2('F', M, M, RWORK(IRVT).asMatrix(M), M,
            WORK(IVT).asMatrix(LDWKVT), LDWKVT);
        zunmbr(
            'P',
            'R',
            'C',
            M,
            M,
            M,
            WORK(IL).asMatrix(LDWRKL),
            LDWRKL,
            WORK(ITAUP),
            WORK(IVT).asMatrix(LDWKVT),
            LDWKVT,
            WORK(NWORK),
            LWORK - NWORK + 1,
            IERR);

        // Multiply right singular vectors of L in WORK(IL) by Q
        // in A, storing result in WORK(IL) and copying to A
        // CWorkspace: need   M*M [VT] + M*M [L]
        // CWorkspace: prefer M*M [VT] + M*N [L]
        // RWorkspace: need   0

        for (I = 1; I <= N; I += CHUNK) {
          BLK = min(N - I + 1, CHUNK);
          zgemm('N', 'N', M, BLK, M, Complex.one, WORK(IVT).asMatrix(M), M,
              A(1, I), LDA, Complex.zero, WORK(IL).asMatrix(LDWRKL), LDWRKL);
          zlacpy('F', M, BLK, WORK(IL).asMatrix(LDWRKL), LDWRKL, A(1, I), LDA);
        }
      } else if (WNTQS) {
        // Path 3t (N >> M, JOBZ='S')
        // M right singular vectors to be computed in VT and
        // M left singular vectors to be computed in U

        IL = 1;

        // WORK(IL) is M by M

        LDWRKL = M;
        ITAU = IL + LDWRKL * M;
        NWORK = ITAU + M;

        // Compute A=L*Q
        // CWorkspace: need   M*M [L] + M [tau] + M    [work]
        // CWorkspace: prefer M*M [L] + M [tau] + M*NB [work]
        // RWorkspace: need   0

        zgelqf(M, N, A, LDA, WORK(ITAU), WORK(NWORK), LWORK - NWORK + 1, IERR);

        // Copy L to WORK(IL), zeroing out above it

        zlacpy('L', M, M, A, LDA, WORK(IL).asMatrix(LDWRKL), LDWRKL);
        zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero,
            WORK(IL + LDWRKL).asMatrix(LDWRKL), LDWRKL);

        // Generate Q in A
        // CWorkspace: need   M*M [L] + M [tau] + M    [work]
        // CWorkspace: prefer M*M [L] + M [tau] + M*NB [work]
        // RWorkspace: need   0

        zunglq(
            M, N, M, A, LDA, WORK(ITAU), WORK(NWORK), LWORK - NWORK + 1, IERR);
        IE = 1;
        ITAUQ = ITAU;
        ITAUP = ITAUQ + M;
        NWORK = ITAUP + M;

        // Bidiagonalize L in WORK(IL)
        // CWorkspace: need   M*M [L] + 2*M [tauq, taup] + M      [work]
        // CWorkspace: prefer M*M [L] + 2*M [tauq, taup] + 2*M*NB [work]
        // RWorkspace: need   M [e]

        zgebrd(M, M, WORK(IL).asMatrix(LDWRKL), LDWRKL, S, RWORK(IE),
            WORK(ITAUQ), WORK(ITAUP), WORK(NWORK), LWORK - NWORK + 1, IERR);

        // Perform bidiagonal SVD, computing left singular vectors
        // of bidiagonal matrix in RWORK(IRU) and computing right
        // singular vectors of bidiagonal matrix in RWORK(IRVT)
        // CWorkspace: need   0
        // RWorkspace: need   M [e] + M*M [RU] + M*M [RVT] + BDSPAC

        IRU = IE + M;
        IRVT = IRU + M * M;
        NRWORK = IRVT + M * M;
        dbdsdc('U', 'I', M, S, RWORK(IE), RWORK(IRU).asMatrix(M), M,
            RWORK(IRVT).asMatrix(M), M, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);

        // Copy real matrix RWORK(IRU) to complex matrix U
        // Overwrite U by left singular vectors of L
        // CWorkspace: need   M*M [L] + 2*M [tauq, taup] + M    [work]
        // CWorkspace: prefer M*M [L] + 2*M [tauq, taup] + M*NB [work]
        // RWorkspace: need   0

        zlacp2('F', M, M, RWORK(IRU).asMatrix(M), M, U, LDU);
        zunmbr('Q', 'L', 'N', M, M, M, WORK(IL).asMatrix(LDWRKL), LDWRKL,
            WORK(ITAUQ), U, LDU, WORK(NWORK), LWORK - NWORK + 1, IERR);

        // Copy real matrix RWORK(IRVT) to complex matrix VT
        // Overwrite VT by left singular vectors of L
        // CWorkspace: need   M*M [L] + 2*M [tauq, taup] + M    [work]
        // CWorkspace: prefer M*M [L] + 2*M [tauq, taup] + M*NB [work]
        // RWorkspace: need   0

        zlacp2('F', M, M, RWORK(IRVT).asMatrix(M), M, VT, LDVT);
        zunmbr('P', 'R', 'C', M, M, M, WORK(IL).asMatrix(LDWRKL), LDWRKL,
            WORK(ITAUP), VT, LDVT, WORK(NWORK), LWORK - NWORK + 1, IERR);

        // Copy VT to WORK(IL), multiply right singular vectors of L
        // in WORK(IL) by Q in A, storing result in VT
        // CWorkspace: need   M*M [L]
        // RWorkspace: need   0

        zlacpy('F', M, M, VT, LDVT, WORK(IL).asMatrix(LDWRKL), LDWRKL);
        zgemm('N', 'N', M, N, M, Complex.one, WORK(IL).asMatrix(LDWRKL), LDWRKL,
            A, LDA, Complex.zero, VT, LDVT);
      } else if (WNTQA) {
        // Path 4t (N >> M, JOBZ='A')
        // N right singular vectors to be computed in VT and
        // M left singular vectors to be computed in U

        IVT = 1;

        // WORK(IVT) is M by M

        LDWKVT = M;
        ITAU = IVT + LDWKVT * M;
        NWORK = ITAU + M;

        // Compute A=L*Q, copying result to VT
        // CWorkspace: need   M*M [VT] + M [tau] + M    [work]
        // CWorkspace: prefer M*M [VT] + M [tau] + M*NB [work]
        // RWorkspace: need   0

        zgelqf(M, N, A, LDA, WORK(ITAU), WORK(NWORK), LWORK - NWORK + 1, IERR);
        zlacpy('U', M, N, A, LDA, VT, LDVT);

        // Generate Q in VT
        // CWorkspace: need   M*M [VT] + M [tau] + N    [work]
        // CWorkspace: prefer M*M [VT] + M [tau] + N*NB [work]
        // RWorkspace: need   0

        zunglq(N, N, M, VT, LDVT, WORK(ITAU), WORK(NWORK), LWORK - NWORK + 1,
            IERR);

        // Produce L in A, zeroing out above it

        zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero, A(1, 2), LDA);
        IE = 1;
        ITAUQ = ITAU;
        ITAUP = ITAUQ + M;
        NWORK = ITAUP + M;

        // Bidiagonalize L in A
        // CWorkspace: need   M*M [VT] + 2*M [tauq, taup] + M      [work]
        // CWorkspace: prefer M*M [VT] + 2*M [tauq, taup] + 2*M*NB [work]
        // RWorkspace: need   M [e]

        zgebrd(M, M, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
            WORK(NWORK), LWORK - NWORK + 1, IERR);

        // Perform bidiagonal SVD, computing left singular vectors
        // of bidiagonal matrix in RWORK(IRU) and computing right
        // singular vectors of bidiagonal matrix in RWORK(IRVT)
        // CWorkspace: need   0
        // RWorkspace: need   M [e] + M*M [RU] + M*M [RVT] + BDSPAC

        IRU = IE + M;
        IRVT = IRU + M * M;
        NRWORK = IRVT + M * M;
        dbdsdc('U', 'I', M, S, RWORK(IE), RWORK(IRU).asMatrix(M), M,
            RWORK(IRVT).asMatrix(M), M, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);

        // Copy real matrix RWORK(IRU) to complex matrix U
        // Overwrite U by left singular vectors of L
        // CWorkspace: need   M*M [VT] + 2*M [tauq, taup] + M    [work]
        // CWorkspace: prefer M*M [VT] + 2*M [tauq, taup] + M*NB [work]
        // RWorkspace: need   0

        zlacp2('F', M, M, RWORK(IRU).asMatrix(M), M, U, LDU);
        zunmbr('Q', 'L', 'N', M, M, M, A, LDA, WORK(ITAUQ), U, LDU, WORK(NWORK),
            LWORK - NWORK + 1, IERR);

        // Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT)
        // Overwrite WORK(IVT) by right singular vectors of L
        // CWorkspace: need   M*M [VT] + 2*M [tauq, taup] + M    [work]
        // CWorkspace: prefer M*M [VT] + 2*M [tauq, taup] + M*NB [work]
        // RWorkspace: need   0

        zlacp2('F', M, M, RWORK(IRVT).asMatrix(M), M,
            WORK(IVT).asMatrix(LDWKVT), LDWKVT);
        zunmbr(
            'P',
            'R',
            'C',
            M,
            M,
            M,
            A,
            LDA,
            WORK(ITAUP),
            WORK(IVT).asMatrix(LDWKVT),
            LDWKVT,
            WORK(NWORK),
            LWORK - NWORK + 1,
            IERR);

        // Multiply right singular vectors of L in WORK(IVT) by
        // Q in VT, storing result in A
        // CWorkspace: need   M*M [VT]
        // RWorkspace: need   0

        zgemm('N', 'N', M, N, M, Complex.one, WORK(IVT).asMatrix(LDWKVT),
            LDWKVT, VT, LDVT, Complex.zero, A, LDA);

        // Copy right singular vectors of A from A to VT

        zlacpy('F', M, N, A, LDA, VT, LDVT);
      }
    } else if (N >= MNTHR2) {
      // MNTHR2 <= N < MNTHR1

      // Path 5t (N >> M, but not as much as MNTHR1)
      // Reduce to bidiagonal form without QR decomposition, use
      // ZUNGBR and matrix multiplication to compute singular vectors

      IE = 1;
      NRWORK = IE + M;
      ITAUQ = 1;
      ITAUP = ITAUQ + M;
      NWORK = ITAUP + M;

      // Bidiagonalize A
      // CWorkspace: need   2*M [tauq, taup] + N        [work]
      // CWorkspace: prefer 2*M [tauq, taup] + (M+N)*NB [work]
      // RWorkspace: need   M [e]

      zgebrd(M, N, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP), WORK(NWORK),
          LWORK - NWORK + 1, IERR);

      if (WNTQN) {
        // Path 5tn (N >> M, JOBZ='N')
        // Compute singular values only
        // CWorkspace: need   0
        // RWorkspace: need   M [e] + BDSPAC

        dbdsdc('L', 'N', M, S, RWORK(IE), DUM.asMatrix(1), 1, DUM.asMatrix(1),
            1, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);
      } else if (WNTQO) {
        IRVT = NRWORK;
        IRU = IRVT + M * M;
        NRWORK = IRU + M * M;
        IVT = NWORK;

        // Path 5to (N >> M, JOBZ='O')
        // Copy A to U, generate Q
        // CWorkspace: need   2*M [tauq, taup] + M    [work]
        // CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
        // RWorkspace: need   0

        zlacpy('L', M, M, A, LDA, U, LDU);
        zungbr('Q', M, M, N, U, LDU, WORK(ITAUQ), WORK(NWORK),
            LWORK - NWORK + 1, IERR);

        // Generate P**H in A
        // CWorkspace: need   2*M [tauq, taup] + M    [work]
        // CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
        // RWorkspace: need   0

        zungbr('P', M, N, M, A, LDA, WORK(ITAUP), WORK(NWORK),
            LWORK - NWORK + 1, IERR);

        LDWKVT = M;
        if (LWORK >= M * N + 3 * M) {
          // WORK( IVT ) is M by N

          NWORK = IVT + LDWKVT * N;
          CHUNK = N;
        } else {
          // WORK( IVT ) is M by CHUNK

          CHUNK = (LWORK - 3 * M) ~/ M;
          NWORK = IVT + LDWKVT * CHUNK;
        }

        // Perform bidiagonal SVD, computing left singular vectors
        // of bidiagonal matrix in RWORK(IRU) and computing right
        // singular vectors of bidiagonal matrix in RWORK(IRVT)
        // CWorkspace: need   0
        // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC

        dbdsdc('L', 'I', M, S, RWORK(IE), RWORK(IRU).asMatrix(M), M,
            RWORK(IRVT).asMatrix(M), M, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);

        // Multiply Q in U by real matrix RWORK(IRVT)
        // storing the result in WORK(IVT), copying to U
        // CWorkspace: need   2*M [tauq, taup] + M*M [VT]
        // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + 2*M*M [rwork]

        zlacrm(M, M, U, LDU, RWORK(IRU).asMatrix(M), M,
            WORK(IVT).asMatrix(LDWKVT), LDWKVT, RWORK(NRWORK));
        zlacpy('F', M, M, WORK(IVT).asMatrix(LDWKVT), LDWKVT, U, LDU);

        // Multiply RWORK(IRVT) by P**H in A, storing the
        // result in WORK(IVT), copying to A
        // CWorkspace: need   2*M [tauq, taup] + M*M [VT]
        // CWorkspace: prefer 2*M [tauq, taup] + M*N [VT]
        // RWorkspace: need   M [e] + M*M [RVT] + 2*M*M [rwork]
        // RWorkspace: prefer M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here

        NRWORK = IRU;
        for (I = 1; I <= N; I += CHUNK) {
          BLK = min(N - I + 1, CHUNK);
          zlarcm(M, BLK, RWORK(IRVT).asMatrix(M), M, A(1, I), LDA,
              WORK(IVT).asMatrix(LDWKVT), LDWKVT, RWORK(NRWORK));
          zlacpy('F', M, BLK, WORK(IVT).asMatrix(LDWKVT), LDWKVT, A(1, I), LDA);
        }
      } else if (WNTQS) {
        // Path 5ts (N >> M, JOBZ='S')
        // Copy A to U, generate Q
        // CWorkspace: need   2*M [tauq, taup] + M    [work]
        // CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
        // RWorkspace: need   0

        zlacpy('L', M, M, A, LDA, U, LDU);
        zungbr('Q', M, M, N, U, LDU, WORK(ITAUQ), WORK(NWORK),
            LWORK - NWORK + 1, IERR);

        // Copy A to VT, generate P**H
        // CWorkspace: need   2*M [tauq, taup] + M    [work]
        // CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
        // RWorkspace: need   0

        zlacpy('U', M, N, A, LDA, VT, LDVT);
        zungbr('P', M, N, M, VT, LDVT, WORK(ITAUP), WORK(NWORK),
            LWORK - NWORK + 1, IERR);

        // Perform bidiagonal SVD, computing left singular vectors
        // of bidiagonal matrix in RWORK(IRU) and computing right
        // singular vectors of bidiagonal matrix in RWORK(IRVT)
        // CWorkspace: need   0
        // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC

        IRVT = NRWORK;
        IRU = IRVT + M * M;
        NRWORK = IRU + M * M;
        dbdsdc('L', 'I', M, S, RWORK(IE), RWORK(IRU).asMatrix(M), M,
            RWORK(IRVT).asMatrix(M), M, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);

        // Multiply Q in U by real matrix RWORK(IRU), storing the
        // result in A, copying to U
        // CWorkspace: need   0
        // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + 2*M*M [rwork]

        zlacrm(M, M, U, LDU, RWORK(IRU).asMatrix(M), M, A, LDA, RWORK(NRWORK));
        zlacpy('F', M, M, A, LDA, U, LDU);

        // Multiply real matrix RWORK(IRVT) by P**H in VT,
        // storing the result in A, copying to VT
        // CWorkspace: need   0
        // RWorkspace: need   M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here

        NRWORK = IRU;
        zlarcm(
            M, N, RWORK(IRVT).asMatrix(M), M, VT, LDVT, A, LDA, RWORK(NRWORK));
        zlacpy('F', M, N, A, LDA, VT, LDVT);
      } else {
        // Path 5ta (N >> M, JOBZ='A')
        // Copy A to U, generate Q
        // CWorkspace: need   2*M [tauq, taup] + M    [work]
        // CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
        // RWorkspace: need   0

        zlacpy('L', M, M, A, LDA, U, LDU);
        zungbr('Q', M, M, N, U, LDU, WORK(ITAUQ), WORK(NWORK),
            LWORK - NWORK + 1, IERR);

        // Copy A to VT, generate P**H
        // CWorkspace: need   2*M [tauq, taup] + N    [work]
        // CWorkspace: prefer 2*M [tauq, taup] + N*NB [work]
        // RWorkspace: need   0

        zlacpy('U', M, N, A, LDA, VT, LDVT);
        zungbr('P', N, N, M, VT, LDVT, WORK(ITAUP), WORK(NWORK),
            LWORK - NWORK + 1, IERR);

        // Perform bidiagonal SVD, computing left singular vectors
        // of bidiagonal matrix in RWORK(IRU) and computing right
        // singular vectors of bidiagonal matrix in RWORK(IRVT)
        // CWorkspace: need   0
        // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC

        IRVT = NRWORK;
        IRU = IRVT + M * M;
        NRWORK = IRU + M * M;
        dbdsdc('L', 'I', M, S, RWORK(IE), RWORK(IRU).asMatrix(M), M,
            RWORK(IRVT).asMatrix(M), M, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);

        // Multiply Q in U by real matrix RWORK(IRU), storing the
        // result in A, copying to U
        // CWorkspace: need   0
        // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + 2*M*M [rwork]

        zlacrm(M, M, U, LDU, RWORK(IRU).asMatrix(M), M, A, LDA, RWORK(NRWORK));
        zlacpy('F', M, M, A, LDA, U, LDU);

        // Multiply real matrix RWORK(IRVT) by P**H in VT,
        // storing the result in A, copying to VT
        // CWorkspace: need   0
        // RWorkspace: need   M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here

        NRWORK = IRU;
        zlarcm(
            M, N, RWORK(IRVT).asMatrix(M), M, VT, LDVT, A, LDA, RWORK(NRWORK));
        zlacpy('F', M, N, A, LDA, VT, LDVT);
      }
    } else {
      // N < MNTHR2

      // Path 6t (N > M, but not much larger)
      // Reduce to bidiagonal form without LQ decomposition
      // Use ZUNMBR to compute singular vectors

      IE = 1;
      NRWORK = IE + M;
      ITAUQ = 1;
      ITAUP = ITAUQ + M;
      NWORK = ITAUP + M;

      // Bidiagonalize A
      // CWorkspace: need   2*M [tauq, taup] + N        [work]
      // CWorkspace: prefer 2*M [tauq, taup] + (M+N)*NB [work]
      // RWorkspace: need   M [e]

      zgebrd(M, N, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP), WORK(NWORK),
          LWORK - NWORK + 1, IERR);
      if (WNTQN) {
        // Path 6tn (N > M, JOBZ='N')
        // Compute singular values only
        // CWorkspace: need   0
        // RWorkspace: need   M [e] + BDSPAC

        dbdsdc('L', 'N', M, S, RWORK(IE), DUM.asMatrix(1), 1, DUM.asMatrix(1),
            1, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);
      } else if (WNTQO) {
        // Path 6to (N > M, JOBZ='O')
        LDWKVT = M;
        IVT = NWORK;
        if (LWORK >= M * N + 3 * M) {
          // WORK( IVT ) is M by N

          zlaset('F', M, N, Complex.zero, Complex.zero,
              WORK(IVT).asMatrix(LDWKVT), LDWKVT);
          NWORK = IVT + LDWKVT * N;
        } else {
          // WORK( IVT ) is M by CHUNK

          CHUNK = (LWORK - 3 * M) ~/ M;
          NWORK = IVT + LDWKVT * CHUNK;
        }

        // Perform bidiagonal SVD, computing left singular vectors
        // of bidiagonal matrix in RWORK(IRU) and computing right
        // singular vectors of bidiagonal matrix in RWORK(IRVT)
        // CWorkspace: need   0
        // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC

        IRVT = NRWORK;
        IRU = IRVT + M * M;
        NRWORK = IRU + M * M;
        dbdsdc('L', 'I', M, S, RWORK(IE), RWORK(IRU).asMatrix(M), M,
            RWORK(IRVT).asMatrix(M), M, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);

        // Copy real matrix RWORK(IRU) to complex matrix U
        // Overwrite U by left singular vectors of A
        // CWorkspace: need   2*M [tauq, taup] + M*M [VT] + M    [work]
        // CWorkspace: prefer 2*M [tauq, taup] + M*M [VT] + M*NB [work]
        // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU]

        zlacp2('F', M, M, RWORK(IRU).asMatrix(M), M, U, LDU);
        zunmbr('Q', 'L', 'N', M, M, N, A, LDA, WORK(ITAUQ), U, LDU, WORK(NWORK),
            LWORK - NWORK + 1, IERR);

        if (LWORK >= M * N + 3 * M) {
          // Path 6to-fast
          // Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT)
          // Overwrite WORK(IVT) by right singular vectors of A,
          // copying to A
          // CWorkspace: need   2*M [tauq, taup] + M*N [VT] + M    [work]
          // CWorkspace: prefer 2*M [tauq, taup] + M*N [VT] + M*NB [work]
          // RWorkspace: need   M [e] + M*M [RVT]

          zlacp2('F', M, M, RWORK(IRVT).asMatrix(M), M,
              WORK(IVT).asMatrix(LDWKVT), LDWKVT);
          zunmbr(
              'P',
              'R',
              'C',
              M,
              N,
              M,
              A,
              LDA,
              WORK(ITAUP),
              WORK(IVT).asMatrix(LDWKVT),
              LDWKVT,
              WORK(NWORK),
              LWORK - NWORK + 1,
              IERR);
          zlacpy('F', M, N, WORK(IVT).asMatrix(LDWKVT), LDWKVT, A, LDA);
        } else {
          // Path 6to-slow
          // Generate P**H in A
          // CWorkspace: need   2*M [tauq, taup] + M*M [VT] + M    [work]
          // CWorkspace: prefer 2*M [tauq, taup] + M*M [VT] + M*NB [work]
          // RWorkspace: need   0

          zungbr('P', M, N, M, A, LDA, WORK(ITAUP), WORK(NWORK),
              LWORK - NWORK + 1, IERR);

          // Multiply Q in A by real matrix RWORK(IRU), storing the
          // result in WORK(IU), copying to A
          // CWorkspace: need   2*M [tauq, taup] + M*M [VT]
          // CWorkspace: prefer 2*M [tauq, taup] + M*N [VT]
          // RWorkspace: need   M [e] + M*M [RVT] + 2*M*M [rwork]
          // RWorkspace: prefer M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here

          NRWORK = IRU;
          for (I = 1; I <= N; I += CHUNK) {
            BLK = min(N - I + 1, CHUNK);
            zlarcm(M, BLK, RWORK(IRVT).asMatrix(M), M, A(1, I), LDA,
                WORK(IVT).asMatrix(LDWKVT), LDWKVT, RWORK(NRWORK));
            zlacpy(
                'F', M, BLK, WORK(IVT).asMatrix(LDWKVT), LDWKVT, A(1, I), LDA);
          }
        }
      } else if (WNTQS) {
        // Path 6ts (N > M, JOBZ='S')
        // Perform bidiagonal SVD, computing left singular vectors
        // of bidiagonal matrix in RWORK(IRU) and computing right
        // singular vectors of bidiagonal matrix in RWORK(IRVT)
        // CWorkspace: need   0
        // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC

        IRVT = NRWORK;
        IRU = IRVT + M * M;
        NRWORK = IRU + M * M;
        dbdsdc('L', 'I', M, S, RWORK(IE), RWORK(IRU).asMatrix(M), M,
            RWORK(IRVT).asMatrix(M), M, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);

        // Copy real matrix RWORK(IRU) to complex matrix U
        // Overwrite U by left singular vectors of A
        // CWorkspace: need   2*M [tauq, taup] + M    [work]
        // CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
        // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU]

        zlacp2('F', M, M, RWORK(IRU).asMatrix(M), M, U, LDU);
        zunmbr('Q', 'L', 'N', M, M, N, A, LDA, WORK(ITAUQ), U, LDU, WORK(NWORK),
            LWORK - NWORK + 1, IERR);

        // Copy real matrix RWORK(IRVT) to complex matrix VT
        // Overwrite VT by right singular vectors of A
        // CWorkspace: need   2*M [tauq, taup] + M    [work]
        // CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
        // RWorkspace: need   M [e] + M*M [RVT]

        zlaset('F', M, N, Complex.zero, Complex.zero, VT, LDVT);
        zlacp2('F', M, M, RWORK(IRVT).asMatrix(M), M, VT, LDVT);
        zunmbr('P', 'R', 'C', M, N, M, A, LDA, WORK(ITAUP), VT, LDVT,
            WORK(NWORK), LWORK - NWORK + 1, IERR);
      } else {
        // Path 6ta (N > M, JOBZ='A')
        // Perform bidiagonal SVD, computing left singular vectors
        // of bidiagonal matrix in RWORK(IRU) and computing right
        // singular vectors of bidiagonal matrix in RWORK(IRVT)
        // CWorkspace: need   0
        // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC

        IRVT = NRWORK;
        IRU = IRVT + M * M;
        NRWORK = IRU + M * M;

        dbdsdc('L', 'I', M, S, RWORK(IE), RWORK(IRU).asMatrix(M), M,
            RWORK(IRVT).asMatrix(M), M, DUM, IDUM, RWORK(NRWORK), IWORK, INFO);

        // Copy real matrix RWORK(IRU) to complex matrix U
        // Overwrite U by left singular vectors of A
        // CWorkspace: need   2*M [tauq, taup] + M    [work]
        // CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
        // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU]

        zlacp2('F', M, M, RWORK(IRU).asMatrix(M), M, U, LDU);
        zunmbr('Q', 'L', 'N', M, M, N, A, LDA, WORK(ITAUQ), U, LDU, WORK(NWORK),
            LWORK - NWORK + 1, IERR);

        // Set all of VT to identity matrix

        zlaset('F', N, N, Complex.zero, Complex.one, VT, LDVT);

        // Copy real matrix RWORK(IRVT) to complex matrix VT
        // Overwrite VT by right singular vectors of A
        // CWorkspace: need   2*M [tauq, taup] + N    [work]
        // CWorkspace: prefer 2*M [tauq, taup] + N*NB [work]
        // RWorkspace: need   M [e] + M*M [RVT]

        zlacp2('F', M, M, RWORK(IRVT).asMatrix(M), M, VT, LDVT);
        zunmbr('P', 'R', 'C', N, N, M, A, LDA, WORK(ITAUP), VT, LDVT,
            WORK(NWORK), LWORK - NWORK + 1, IERR);
      }
    }
  }

  // Undo scaling if necessary

  if (ISCL == 1) {
    if (ANRM > BIGNUM) {
      dlascl('G', 0, 0, BIGNUM, ANRM, MINMN, 1, S.asMatrix(MINMN), MINMN, IERR);
    }
    if (INFO.value != 0 && ANRM > BIGNUM) {
      dlascl('G', 0, 0, BIGNUM, ANRM, MINMN - 1, 1, RWORK(IE).asMatrix(MINMN),
          MINMN, IERR);
    }
    if (ANRM < SMLNUM) {
      dlascl('G', 0, 0, SMLNUM, ANRM, MINMN, 1, S.asMatrix(MINMN), MINMN, IERR);
    }
    if (INFO.value != 0 && ANRM < SMLNUM) {
      dlascl('G', 0, 0, SMLNUM, ANRM, MINMN - 1, 1, RWORK(IE).asMatrix(MINMN),
          MINMN, IERR);
    }
  }

  // Return optimal workspace in WORK(1)

  WORK[1] = droundup_lwork(MAXWRK).toComplex();
}
