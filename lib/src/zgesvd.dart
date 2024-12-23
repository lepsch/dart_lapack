// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/dlascl.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zgeqrf.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zbdsqr.dart';
import 'package:dart_lapack/src/zgebrd.dart';
import 'package:dart_lapack/src/zgelqf.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlascl.dart';
import 'package:dart_lapack/src/zlaset.dart';
import 'package:dart_lapack/src/zungbr.dart';
import 'package:dart_lapack/src/zunglq.dart';
import 'package:dart_lapack/src/zungqr.dart';
import 'package:dart_lapack/src/zunmbr.dart';

void zgesvd(
  final String JOBU,
  final String JOBVT,
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
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final U = U_.having(ld: LDU);
  final VT = VT_.having(ld: LDVT);
  final S = S_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();

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
  int BLK,
      CHUNK,
      I,
      IE = 0,
      IR,
      IRWORK,
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
  int LWORK_ZGEQRF,
      LWORK_ZUNGQR_N,
      LWORK_ZUNGQR_M,
      LWORK_ZGEBRD,
      LWORK_ZUNGBR_P,
      LWORK_ZUNGBR_Q,
      LWORK_ZGELQF,
      LWORK_ZUNGLQ_N,
      LWORK_ZUNGLQ_M;
  double ANRM, BIGNUM, EPS, SMLNUM;
  final DUM = Array<double>(1);
  final CDUM = Array<Complex>(1);
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
  //   CWorkspace refers to complex workspace, and RWorkspace to
  //   real workspace. NB refers to the optimal block size for the
  //   immediately following subroutine, as returned by ILAENV.)

  if (INFO.value == 0) {
    MINWRK = 1;
    MAXWRK = 1;
    if (M >= N && MINMN > 0) {
      // Space needed for ZBDSQR is BDSPAC = 5*N

      MNTHR = ilaenv(6, 'ZGESVD', JOBU + JOBVT, M, N, 0, 0);
      // Compute space needed for ZGEQRF
      zgeqrf(M, N, A, LDA, CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZGEQRF = CDUM[1].toInt();
      // Compute space needed for ZUNGQR
      zungqr(M, N, N, A, LDA, CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZUNGQR_N = CDUM[1].toInt();
      zungqr(M, M, N, A, LDA, CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZUNGQR_M = CDUM[1].toInt();
      // Compute space needed for ZGEBRD
      zgebrd(N, N, A, LDA, S, DUM(1), CDUM(1), CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZGEBRD = CDUM[1].toInt();
      // Compute space needed for ZUNGBR
      zungbr('P', N, N, N, A, LDA, CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZUNGBR_P = CDUM[1].toInt();
      zungbr('Q', N, N, N, A, LDA, CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZUNGBR_Q = CDUM[1].toInt();

      if (M >= MNTHR) {
        if (WNTUN) {
          // Path 1 (M much larger than N, JOBU='N')

          MAXWRK = N + LWORK_ZGEQRF;
          MAXWRK = max(MAXWRK, 2 * N + LWORK_ZGEBRD);
          if (WNTVO || WNTVAS) MAXWRK = max(MAXWRK, 2 * N + LWORK_ZUNGBR_P);
          MINWRK = 3 * N;
        } else if (WNTUO && WNTVN) {
          // Path 2 (M much larger than N, JOBU='O', JOBVT='N')

          WRKBL = N + LWORK_ZGEQRF;
          WRKBL = max(WRKBL, N + LWORK_ZUNGQR_N);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZGEBRD);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZUNGBR_Q);
          MAXWRK = max(N * N + WRKBL, N * N + M * N);
          MINWRK = 2 * N + M;
        } else if (WNTUO && WNTVAS) {
          // Path 3 (M much larger than N, JOBU='O', JOBVT='S' or
          // 'A')

          WRKBL = N + LWORK_ZGEQRF;
          WRKBL = max(WRKBL, N + LWORK_ZUNGQR_N);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZGEBRD);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZUNGBR_Q);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZUNGBR_P);
          MAXWRK = max(N * N + WRKBL, N * N + M * N);
          MINWRK = 2 * N + M;
        } else if (WNTUS && WNTVN) {
          // Path 4 (M much larger than N, JOBU='S', JOBVT='N')

          WRKBL = N + LWORK_ZGEQRF;
          WRKBL = max(WRKBL, N + LWORK_ZUNGQR_N);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZGEBRD);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZUNGBR_Q);
          MAXWRK = N * N + WRKBL;
          MINWRK = 2 * N + M;
        } else if (WNTUS && WNTVO) {
          // Path 5 (M much larger than N, JOBU='S', JOBVT='O')

          WRKBL = N + LWORK_ZGEQRF;
          WRKBL = max(WRKBL, N + LWORK_ZUNGQR_N);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZGEBRD);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZUNGBR_Q);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZUNGBR_P);
          MAXWRK = 2 * N * N + WRKBL;
          MINWRK = 2 * N + M;
        } else if (WNTUS && WNTVAS) {
          // Path 6 (M much larger than N, JOBU='S', JOBVT='S' or
          // 'A')

          WRKBL = N + LWORK_ZGEQRF;
          WRKBL = max(WRKBL, N + LWORK_ZUNGQR_N);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZGEBRD);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZUNGBR_Q);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZUNGBR_P);
          MAXWRK = N * N + WRKBL;
          MINWRK = 2 * N + M;
        } else if (WNTUA && WNTVN) {
          // Path 7 (M much larger than N, JOBU='A', JOBVT='N')

          WRKBL = N + LWORK_ZGEQRF;
          WRKBL = max(WRKBL, N + LWORK_ZUNGQR_M);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZGEBRD);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZUNGBR_Q);
          MAXWRK = N * N + WRKBL;
          MINWRK = 2 * N + M;
        } else if (WNTUA && WNTVO) {
          // Path 8 (M much larger than N, JOBU='A', JOBVT='O')

          WRKBL = N + LWORK_ZGEQRF;
          WRKBL = max(WRKBL, N + LWORK_ZUNGQR_M);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZGEBRD);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZUNGBR_Q);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZUNGBR_P);
          MAXWRK = 2 * N * N + WRKBL;
          MINWRK = 2 * N + M;
        } else if (WNTUA && WNTVAS) {
          // Path 9 (M much larger than N, JOBU='A', JOBVT='S' or
          // 'A')

          WRKBL = N + LWORK_ZGEQRF;
          WRKBL = max(WRKBL, N + LWORK_ZUNGQR_M);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZGEBRD);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZUNGBR_Q);
          WRKBL = max(WRKBL, 2 * N + LWORK_ZUNGBR_P);
          MAXWRK = N * N + WRKBL;
          MINWRK = 2 * N + M;
        }
      } else {
        // Path 10 (M at least N, but not much larger)

        zgebrd(M, N, A, LDA, S, DUM(1), CDUM(1), CDUM(1), CDUM(1), -1, IERR);
        LWORK_ZGEBRD = CDUM[1].toInt();
        MAXWRK = 2 * N + LWORK_ZGEBRD;
        if (WNTUS || WNTUO) {
          zungbr('Q', M, N, N, A, LDA, CDUM(1), CDUM(1), -1, IERR);
          LWORK_ZUNGBR_Q = CDUM[1].toInt();
          MAXWRK = max(MAXWRK, 2 * N + LWORK_ZUNGBR_Q);
        }
        if (WNTUA) {
          zungbr('Q', M, M, N, A, LDA, CDUM(1), CDUM(1), -1, IERR);
          LWORK_ZUNGBR_Q = CDUM[1].toInt();
          MAXWRK = max(MAXWRK, 2 * N + LWORK_ZUNGBR_Q);
        }
        if (!WNTVN) {
          MAXWRK = max(MAXWRK, 2 * N + LWORK_ZUNGBR_P);
        }
        MINWRK = 2 * N + M;
      }
    } else if (MINMN > 0) {
      // Space needed for ZBDSQR is BDSPAC = 5*M

      MNTHR = ilaenv(6, 'ZGESVD', JOBU + JOBVT, M, N, 0, 0);
      // Compute space needed for ZGELQF
      zgelqf(M, N, A, LDA, CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZGELQF = CDUM[1].toInt();
      // Compute space needed for ZUNGLQ
      zunglq(N, N, M, CDUM(1).asMatrix(N), N, CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZUNGLQ_N = CDUM[1].toInt();
      zunglq(M, N, M, A, LDA, CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZUNGLQ_M = CDUM[1].toInt();
      // Compute space needed for ZGEBRD
      zgebrd(M, M, A, LDA, S, DUM(1), CDUM(1), CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZGEBRD = CDUM[1].toInt();
      // Compute space needed for ZUNGBR P
      zungbr('P', M, M, M, A, N, CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZUNGBR_P = CDUM[1].toInt();
      // Compute space needed for ZUNGBR Q
      zungbr('Q', M, M, M, A, N, CDUM(1), CDUM(1), -1, IERR);
      LWORK_ZUNGBR_Q = CDUM[1].toInt();
      if (N >= MNTHR) {
        if (WNTVN) {
          // Path 1t(N much larger than M, JOBVT='N')

          MAXWRK = M + LWORK_ZGELQF;
          MAXWRK = max(MAXWRK, 2 * M + LWORK_ZGEBRD);
          if (WNTUO || WNTUAS) MAXWRK = max(MAXWRK, 2 * M + LWORK_ZUNGBR_Q);
          MINWRK = 3 * M;
        } else if (WNTVO && WNTUN) {
          // Path 2t(N much larger than M, JOBU='N', JOBVT='O')

          WRKBL = M + LWORK_ZGELQF;
          WRKBL = max(WRKBL, M + LWORK_ZUNGLQ_M);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZGEBRD);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZUNGBR_P);
          MAXWRK = max(M * M + WRKBL, M * M + M * N);
          MINWRK = 2 * M + N;
        } else if (WNTVO && WNTUAS) {
          // Path 3t(N much larger than M, JOBU='S' or 'A',
          // JOBVT='O')

          WRKBL = M + LWORK_ZGELQF;
          WRKBL = max(WRKBL, M + LWORK_ZUNGLQ_M);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZGEBRD);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZUNGBR_P);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZUNGBR_Q);
          MAXWRK = max(M * M + WRKBL, M * M + M * N);
          MINWRK = 2 * M + N;
        } else if (WNTVS && WNTUN) {
          // Path 4t(N much larger than M, JOBU='N', JOBVT='S')

          WRKBL = M + LWORK_ZGELQF;
          WRKBL = max(WRKBL, M + LWORK_ZUNGLQ_M);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZGEBRD);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZUNGBR_P);
          MAXWRK = M * M + WRKBL;
          MINWRK = 2 * M + N;
        } else if (WNTVS && WNTUO) {
          // Path 5t(N much larger than M, JOBU='O', JOBVT='S')

          WRKBL = M + LWORK_ZGELQF;
          WRKBL = max(WRKBL, M + LWORK_ZUNGLQ_M);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZGEBRD);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZUNGBR_P);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZUNGBR_Q);
          MAXWRK = 2 * M * M + WRKBL;
          MINWRK = 2 * M + N;
        } else if (WNTVS && WNTUAS) {
          // Path 6t(N much larger than M, JOBU='S' or 'A',
          // JOBVT='S')

          WRKBL = M + LWORK_ZGELQF;
          WRKBL = max(WRKBL, M + LWORK_ZUNGLQ_M);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZGEBRD);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZUNGBR_P);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZUNGBR_Q);
          MAXWRK = M * M + WRKBL;
          MINWRK = 2 * M + N;
        } else if (WNTVA && WNTUN) {
          // Path 7t(N much larger than M, JOBU='N', JOBVT='A')

          WRKBL = M + LWORK_ZGELQF;
          WRKBL = max(WRKBL, M + LWORK_ZUNGLQ_N);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZGEBRD);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZUNGBR_P);
          MAXWRK = M * M + WRKBL;
          MINWRK = 2 * M + N;
        } else if (WNTVA && WNTUO) {
          // Path 8t(N much larger than M, JOBU='O', JOBVT='A')

          WRKBL = M + LWORK_ZGELQF;
          WRKBL = max(WRKBL, M + LWORK_ZUNGLQ_N);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZGEBRD);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZUNGBR_P);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZUNGBR_Q);
          MAXWRK = 2 * M * M + WRKBL;
          MINWRK = 2 * M + N;
        } else if (WNTVA && WNTUAS) {
          // Path 9t(N much larger than M, JOBU='S' or 'A',
          // JOBVT='A')

          WRKBL = M + LWORK_ZGELQF;
          WRKBL = max(WRKBL, M + LWORK_ZUNGLQ_N);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZGEBRD);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZUNGBR_P);
          WRKBL = max(WRKBL, 2 * M + LWORK_ZUNGBR_Q);
          MAXWRK = M * M + WRKBL;
          MINWRK = 2 * M + N;
        }
      } else {
        // Path 10t(N greater than M, but not much larger)

        zgebrd(M, N, A, LDA, S, DUM(1), CDUM(1), CDUM(1), CDUM(1), -1, IERR);
        LWORK_ZGEBRD = CDUM[1].toInt();
        MAXWRK = 2 * M + LWORK_ZGEBRD;
        if (WNTVS || WNTVO) {
          // Compute space needed for ZUNGBR P
          zungbr('P', M, N, M, A, N, CDUM(1), CDUM(1), -1, IERR);
          LWORK_ZUNGBR_P = CDUM[1].toInt();
          MAXWRK = max(MAXWRK, 2 * M + LWORK_ZUNGBR_P);
        }
        if (WNTVA) {
          zungbr('P', N, N, M, A, N, CDUM(1), CDUM(1), -1, IERR);
          LWORK_ZUNGBR_P = CDUM[1].toInt();
          MAXWRK = max(MAXWRK, 2 * M + LWORK_ZUNGBR_P);
        }
        if (!WNTUN) {
          MAXWRK = max(MAXWRK, 2 * M + LWORK_ZUNGBR_Q);
        }
        MINWRK = 2 * M + N;
      }
    }
    MAXWRK = max(MAXWRK, MINWRK);
    WORK[1] = MAXWRK.toComplex();

    if (LWORK < MINWRK && !LQUERY) {
      INFO.value = -13;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZGESVD', -INFO.value);
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

    if (M >= MNTHR) {
      if (WNTUN) {
        // Path 1 (M much larger than N, JOBU='N')
        // No left singular vectors to be computed

        ITAU = 1;
        IWORK = ITAU + N;

        // Compute A=Q*R
        // (CWorkspace: need 2*N, prefer N+N*NB)
        // (RWorkspace: need 0)

        zgeqrf(M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

        // Zero out below R

        if (N > 1) {
          zlaset('L', N - 1, N - 1, Complex.zero, Complex.zero, A(2, 1), LDA);
        }
        IE = 1;
        ITAUQ = 1;
        ITAUP = ITAUQ + N;
        IWORK = ITAUP + N;

        // Bidiagonalize R in A
        // (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
        // (RWorkspace: need N)

        zgebrd(N, N, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
            WORK(IWORK), LWORK - IWORK + 1, IERR);
        NCVT = 0;
        if (WNTVO || WNTVAS) {
          // If right singular vectors desired, generate P'.
          // (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
          // (RWorkspace: 0)

          zungbr('P', N, N, N, A, LDA, WORK(ITAUP), WORK(IWORK),
              LWORK - IWORK + 1, IERR);
          NCVT = N;
        }
        IRWORK = IE + N;

        // Perform bidiagonal QR iteration, computing right
        // singular vectors of A in A if desired
        // (CWorkspace: 0)
        // (RWorkspace: need BDSPAC)

        zbdsqr('U', N, NCVT, 0, 0, S, RWORK(IE), A, LDA, CDUM.asMatrix(1), 1,
            CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);

        // If right singular vectors desired in VT, copy them there

        if (WNTVAS) zlacpy('F', N, N, A, LDA, VT, LDVT);
      } else if (WNTUO && WNTVN) {
        // Path 2 (M much larger than N, JOBU='O', JOBVT='N')
        // N left singular vectors to be overwritten on A and
        // no right singular vectors to be computed

        if (LWORK >= N * N + 3 * N) {
          // Sufficient workspace for a fast algorithm

          IR = 1;
          if (LWORK >= max(WRKBL, LDA * N) + LDA * N) {
            // WORK(IU) is LDA by N, WORK(IR) is LDA by N

            LDWRKU = LDA;
            LDWRKR = LDA;
          } else if (LWORK >= max(WRKBL, LDA * N) + N * N) {
            // WORK(IU) is LDA by N, WORK(IR) is N by N

            LDWRKU = LDA;
            LDWRKR = N;
          } else {
            // WORK(IU) is LDWRKU by N, WORK(IR) is N by N

            LDWRKU = (LWORK - N * N) ~/ N;
            LDWRKR = N;
          }
          ITAU = IR + LDWRKR * N;
          IWORK = ITAU + N;

          // Compute A=Q*R
          // (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
          // (RWorkspace: 0)

          zgeqrf(
              M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Copy R to WORK(IR) and zero out below it

          zlacpy('U', N, N, A, LDA, WORK(IR).asMatrix(LDWRKR), LDWRKR);
          zlaset('L', N - 1, N - 1, Complex.zero, Complex.zero,
              WORK(IR + 1).asMatrix(LDWRKR), LDWRKR);

          // Generate Q in A
          // (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
          // (RWorkspace: 0)

          zungqr(M, N, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
              IERR);
          IE = 1;
          ITAUQ = ITAU;
          ITAUP = ITAUQ + N;
          IWORK = ITAUP + N;

          // Bidiagonalize R in WORK(IR)
          // (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
          // (RWorkspace: need N)

          zgebrd(N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, S, RWORK(IE),
              WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Generate left vectors bidiagonalizing R
          // (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
          // (RWorkspace: need 0)

          zungbr('Q', N, N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUQ),
              WORK(IWORK), LWORK - IWORK + 1, IERR);
          IRWORK = IE + N;

          // Perform bidiagonal QR iteration, computing left
          // singular vectors of R in WORK(IR)
          // (CWorkspace: need N*N)
          // (RWorkspace: need BDSPAC)

          zbdsqr(
              'U',
              N,
              0,
              N,
              0,
              S,
              RWORK(IE),
              CDUM.asMatrix(1),
              1,
              WORK(IR).asMatrix(LDWRKR),
              LDWRKR,
              CDUM.asMatrix(1),
              1,
              RWORK(IRWORK),
              INFO);
          IU = ITAUQ;

          // Multiply Q in A by left singular vectors of R in
          // WORK(IR), storing result in WORK(IU) and copying to A
          // (CWorkspace: need N*N+N, prefer N*N+M*N)
          // (RWorkspace: 0)

          for (I = 1; I <= M; I += LDWRKU) {
            CHUNK = min(M - I + 1, LDWRKU);
            zgemm(
                'N',
                'N',
                CHUNK,
                N,
                N,
                Complex.one,
                A(I, 1),
                LDA,
                WORK(IR).asMatrix(LDWRKR),
                LDWRKR,
                Complex.zero,
                WORK(IU).asMatrix(LDWRKU),
                LDWRKU);
            zlacpy(
                'F', CHUNK, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, A(I, 1), LDA);
          }
        } else {
          // Insufficient workspace for a fast algorithm

          IE = 1;
          ITAUQ = 1;
          ITAUP = ITAUQ + N;
          IWORK = ITAUP + N;

          // Bidiagonalize A
          // (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB)
          // (RWorkspace: N)

          zgebrd(M, N, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
              WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Generate left vectors bidiagonalizing A
          // (CWorkspace: need 3*N, prefer 2*N+N*NB)
          // (RWorkspace: 0)

          zungbr('Q', M, N, N, A, LDA, WORK(ITAUQ), WORK(IWORK),
              LWORK - IWORK + 1, IERR);
          IRWORK = IE + N;

          // Perform bidiagonal QR iteration, computing left
          // singular vectors of A in A
          // (CWorkspace: need 0)
          // (RWorkspace: need BDSPAC)

          zbdsqr('U', N, 0, M, 0, S, RWORK(IE), CDUM.asMatrix(1), 1, A, LDA,
              CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
        }
      } else if (WNTUO && WNTVAS) {
        // Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A')
        // N left singular vectors to be overwritten on A and
        // N right singular vectors to be computed in VT

        if (LWORK >= N * N + 3 * N) {
          // Sufficient workspace for a fast algorithm

          IR = 1;
          if (LWORK >= max(WRKBL, LDA * N) + LDA * N) {
            // WORK(IU) is LDA by N and WORK(IR) is LDA by N

            LDWRKU = LDA;
            LDWRKR = LDA;
          } else if (LWORK >= max(WRKBL, LDA * N) + N * N) {
            // WORK(IU) is LDA by N and WORK(IR) is N by N

            LDWRKU = LDA;
            LDWRKR = N;
          } else {
            // WORK(IU) is LDWRKU by N and WORK(IR) is N by N

            LDWRKU = (LWORK - N * N) ~/ N;
            LDWRKR = N;
          }
          ITAU = IR + LDWRKR * N;
          IWORK = ITAU + N;

          // Compute A=Q*R
          // (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
          // (RWorkspace: 0)

          zgeqrf(
              M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Copy R to VT, zeroing out below it

          zlacpy('U', N, N, A, LDA, VT, LDVT);
          if (N > 1) {
            zlaset(
                'L', N - 1, N - 1, Complex.zero, Complex.zero, VT(2, 1), LDVT);
          }

          // Generate Q in A
          // (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
          // (RWorkspace: 0)

          zungqr(M, N, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
              IERR);
          IE = 1;
          ITAUQ = ITAU;
          ITAUP = ITAUQ + N;
          IWORK = ITAUP + N;

          // Bidiagonalize R in VT, copying result to WORK(IR)
          // (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
          // (RWorkspace: need N)

          zgebrd(N, N, VT, LDVT, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
              WORK(IWORK), LWORK - IWORK + 1, IERR);
          zlacpy('L', N, N, VT, LDVT, WORK(IR).asMatrix(LDWRKR), LDWRKR);

          // Generate left vectors bidiagonalizing R in WORK(IR)
          // (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
          // (RWorkspace: 0)

          zungbr('Q', N, N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUQ),
              WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Generate right vectors bidiagonalizing R in VT
          // (CWorkspace: need N*N+3*N-1, prefer N*N+2*N+(N-1)*NB)
          // (RWorkspace: 0)

          zungbr('P', N, N, N, VT, LDVT, WORK(ITAUP), WORK(IWORK),
              LWORK - IWORK + 1, IERR);
          IRWORK = IE + N;

          // Perform bidiagonal QR iteration, computing left
          // singular vectors of R in WORK(IR) and computing right
          // singular vectors of R in VT
          // (CWorkspace: need N*N)
          // (RWorkspace: need BDSPAC)

          zbdsqr(
              'U',
              N,
              N,
              N,
              0,
              S,
              RWORK(IE),
              VT,
              LDVT,
              WORK(IR).asMatrix(LDWRKR),
              LDWRKR,
              CDUM.asMatrix(1),
              1,
              RWORK(IRWORK),
              INFO);
          IU = ITAUQ;

          // Multiply Q in A by left singular vectors of R in
          // WORK(IR), storing result in WORK(IU) and copying to A
          // (CWorkspace: need N*N+N, prefer N*N+M*N)
          // (RWorkspace: 0)

          for (I = 1; I <= M; I += LDWRKU) {
            CHUNK = min(M - I + 1, LDWRKU);
            zgemm(
                'N',
                'N',
                CHUNK,
                N,
                N,
                Complex.one,
                A(I, 1),
                LDA,
                WORK(IR).asMatrix(LDWRKR),
                LDWRKR,
                Complex.zero,
                WORK(IU).asMatrix(LDWRKU),
                LDWRKU);
            zlacpy(
                'F', CHUNK, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, A(I, 1), LDA);
          }
        } else {
          // Insufficient workspace for a fast algorithm

          ITAU = 1;
          IWORK = ITAU + N;

          // Compute A=Q*R
          // (CWorkspace: need 2*N, prefer N+N*NB)
          // (RWorkspace: 0)

          zgeqrf(
              M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Copy R to VT, zeroing out below it

          zlacpy('U', N, N, A, LDA, VT, LDVT);
          if (N > 1) {
            zlaset(
                'L', N - 1, N - 1, Complex.zero, Complex.zero, VT(2, 1), LDVT);
          }

          // Generate Q in A
          // (CWorkspace: need 2*N, prefer N+N*NB)
          // (RWorkspace: 0)

          zungqr(M, N, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
              IERR);
          IE = 1;
          ITAUQ = ITAU;
          ITAUP = ITAUQ + N;
          IWORK = ITAUP + N;

          // Bidiagonalize R in VT
          // (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
          // (RWorkspace: N)

          zgebrd(N, N, VT, LDVT, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
              WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Multiply Q in A by left vectors bidiagonalizing R
          // (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
          // (RWorkspace: 0)

          zunmbr('Q', 'R', 'N', M, N, N, VT, LDVT, WORK(ITAUQ), A, LDA,
              WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Generate right vectors bidiagonalizing R in VT
          // (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
          // (RWorkspace: 0)

          zungbr('P', N, N, N, VT, LDVT, WORK(ITAUP), WORK(IWORK),
              LWORK - IWORK + 1, IERR);
          IRWORK = IE + N;

          // Perform bidiagonal QR iteration, computing left
          // singular vectors of A in A and computing right
          // singular vectors of A in VT
          // (CWorkspace: 0)
          // (RWorkspace: need BDSPAC)

          zbdsqr('U', N, N, M, 0, S, RWORK(IE), VT, LDVT, A, LDA,
              CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
        }
      } else if (WNTUS) {
        if (WNTVN) {
          // Path 4 (M much larger than N, JOBU='S', JOBVT='N')
          // N left singular vectors to be computed in U and
          // no right singular vectors to be computed

          if (LWORK >= N * N + 3 * N) {
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
            // (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
            // (RWorkspace: 0)

            zgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Copy R to WORK(IR), zeroing out below it

            zlacpy('U', N, N, A, LDA, WORK(IR).asMatrix(LDWRKR), LDWRKR);
            zlaset('L', N - 1, N - 1, Complex.zero, Complex.zero,
                WORK(IR + 1).asMatrix(LDWRKR), LDWRKR);

            // Generate Q in A
            // (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
            // (RWorkspace: 0)

            zungqr(M, N, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Bidiagonalize R in WORK(IR)
            // (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
            // (RWorkspace: need N)

            zgebrd(N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, S, RWORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate left vectors bidiagonalizing R in WORK(IR)
            // (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
            // (RWorkspace: 0)

            zungbr('Q', N, N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUQ),
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IRWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of R in WORK(IR)
            // (CWorkspace: need N*N)
            // (RWorkspace: need BDSPAC)

            zbdsqr(
                'U',
                N,
                0,
                N,
                0,
                S,
                RWORK(IE),
                CDUM.asMatrix(1),
                1,
                WORK(IR).asMatrix(LDWRKR),
                LDWRKR,
                CDUM.asMatrix(1),
                1,
                RWORK(IRWORK),
                INFO);

            // Multiply Q in A by left singular vectors of R in
            // WORK(IR), storing result in U
            // (CWorkspace: need N*N)
            // (RWorkspace: 0)

            zgemm('N', 'N', M, N, N, Complex.one, A, LDA,
                WORK(IR).asMatrix(LDWRKR), LDWRKR, Complex.zero, U, LDU);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + N;

            // Compute A=Q*R, copying result to U
            // (CWorkspace: need 2*N, prefer N+N*NB)
            // (RWorkspace: 0)

            zgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('L', M, N, A, LDA, U, LDU);

            // Generate Q in U
            // (CWorkspace: need 2*N, prefer N+N*NB)
            // (RWorkspace: 0)

            zungqr(M, N, N, U, LDU, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Zero out below R in A

            if (N > 1) {
              zlaset(
                  'L', N - 1, N - 1, Complex.zero, Complex.zero, A(2, 1), LDA);
            }

            // Bidiagonalize R in A
            // (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
            // (RWorkspace: need N)

            zgebrd(N, N, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply Q in U by left vectors bidiagonalizing R
            // (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
            // (RWorkspace: 0)

            zunmbr('Q', 'R', 'N', M, N, N, A, LDA, WORK(ITAUQ), U, LDU,
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IRWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of A in U
            // (CWorkspace: 0)
            // (RWorkspace: need BDSPAC)

            zbdsqr('U', N, 0, M, 0, S, RWORK(IE), CDUM.asMatrix(1), 1, U, LDU,
                CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
          }
        } else if (WNTVO) {
          // Path 5 (M much larger than N, JOBU='S', JOBVT='O')
          // N left singular vectors to be computed in U and
          // N right singular vectors to be overwritten on A

          if (LWORK >= 2 * N * N + 3 * N) {
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
            // (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)
            // (RWorkspace: 0)

            zgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Copy R to WORK(IU), zeroing out below it

            zlacpy('U', N, N, A, LDA, WORK(IU).asMatrix(LDWRKU), LDWRKU);
            zlaset('L', N - 1, N - 1, Complex.zero, Complex.zero,
                WORK(IU + 1).asMatrix(LDWRKU), LDWRKU);

            // Generate Q in A
            // (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)
            // (RWorkspace: 0)

            zungqr(M, N, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Bidiagonalize R in WORK(IU), copying result to
            // WORK(IR)
            // (CWorkspace: need   2*N*N+3*N,
            //              prefer 2*N*N+2*N+2*N*NB)
            // (RWorkspace: need   N)

            zgebrd(N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, S, RWORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('U', N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU,
                WORK(IR).asMatrix(LDWRKR), LDWRKR);

            // Generate left bidiagonalizing vectors in WORK(IU)
            // (CWorkspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB)
            // (RWorkspace: 0)

            zungbr('Q', N, N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, WORK(ITAUQ),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate right bidiagonalizing vectors in WORK(IR)
            // (CWorkspace: need   2*N*N+3*N-1,
            //              prefer 2*N*N+2*N+(N-1)*NB)
            // (RWorkspace: 0)

            zungbr('P', N, N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IRWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of R in WORK(IU) and computing
            // right singular vectors of R in WORK(IR)
            // (CWorkspace: need 2*N*N)
            // (RWorkspace: need BDSPAC)

            zbdsqr(
                'U',
                N,
                N,
                N,
                0,
                S,
                RWORK(IE),
                WORK(IR).asMatrix(LDWRKR),
                LDWRKR,
                WORK(IU).asMatrix(LDWRKU),
                LDWRKU,
                CDUM.asMatrix(1),
                1,
                RWORK(IRWORK),
                INFO);

            // Multiply Q in A by left singular vectors of R in
            // WORK(IU), storing result in U
            // (CWorkspace: need N*N)
            // (RWorkspace: 0)

            zgemm('N', 'N', M, N, N, Complex.one, A, LDA,
                WORK(IU).asMatrix(LDWRKU), LDWRKU, Complex.zero, U, LDU);

            // Copy right singular vectors of R to A
            // (CWorkspace: need N*N)
            // (RWorkspace: 0)

            zlacpy('F', N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, A, LDA);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + N;

            // Compute A=Q*R, copying result to U
            // (CWorkspace: need 2*N, prefer N+N*NB)
            // (RWorkspace: 0)

            zgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('L', M, N, A, LDA, U, LDU);

            // Generate Q in U
            // (CWorkspace: need 2*N, prefer N+N*NB)
            // (RWorkspace: 0)

            zungqr(M, N, N, U, LDU, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Zero out below R in A

            if (N > 1) {
              zlaset(
                  'L', N - 1, N - 1, Complex.zero, Complex.zero, A(2, 1), LDA);
            }

            // Bidiagonalize R in A
            // (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
            // (RWorkspace: need N)

            zgebrd(N, N, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply Q in U by left vectors bidiagonalizing R
            // (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
            // (RWorkspace: 0)

            zunmbr('Q', 'R', 'N', M, N, N, A, LDA, WORK(ITAUQ), U, LDU,
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate right vectors bidiagonalizing R in A
            // (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
            // (RWorkspace: 0)

            zungbr('P', N, N, N, A, LDA, WORK(ITAUP), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IRWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of A in U and computing right
            // singular vectors of A in A
            // (CWorkspace: 0)
            // (RWorkspace: need BDSPAC)

            zbdsqr('U', N, N, M, 0, S, RWORK(IE), A, LDA, U, LDU,
                CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
          }
        } else if (WNTVAS) {
          // Path 6 (M much larger than N, JOBU='S', JOBVT='S'
          //         or 'A')
          // N left singular vectors to be computed in U and
          // N right singular vectors to be computed in VT

          if (LWORK >= N * N + 3 * N) {
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
            // (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
            // (RWorkspace: 0)

            zgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Copy R to WORK(IU), zeroing out below it

            zlacpy('U', N, N, A, LDA, WORK(IU).asMatrix(LDWRKU), LDWRKU);
            zlaset('L', N - 1, N - 1, Complex.zero, Complex.zero,
                WORK(IU + 1).asMatrix(LDWRKU), LDWRKU);

            // Generate Q in A
            // (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
            // (RWorkspace: 0)

            zungqr(M, N, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Bidiagonalize R in WORK(IU), copying result to VT
            // (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
            // (RWorkspace: need N)

            zgebrd(N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, S, RWORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('U', N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, VT, LDVT);

            // Generate left bidiagonalizing vectors in WORK(IU)
            // (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
            // (RWorkspace: 0)

            zungbr('Q', N, N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, WORK(ITAUQ),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate right bidiagonalizing vectors in VT
            // (CWorkspace: need   N*N+3*N-1,
            //              prefer N*N+2*N+(N-1)*NB)
            // (RWorkspace: 0)

            zungbr('P', N, N, N, VT, LDVT, WORK(ITAUP), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IRWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of R in WORK(IU) and computing
            // right singular vectors of R in VT
            // (CWorkspace: need N*N)
            // (RWorkspace: need BDSPAC)

            zbdsqr(
                'U',
                N,
                N,
                N,
                0,
                S,
                RWORK(IE),
                VT,
                LDVT,
                WORK(IU).asMatrix(LDWRKU),
                LDWRKU,
                CDUM.asMatrix(1),
                1,
                RWORK(IRWORK),
                INFO);

            // Multiply Q in A by left singular vectors of R in
            // WORK(IU), storing result in U
            // (CWorkspace: need N*N)
            // (RWorkspace: 0)

            zgemm('N', 'N', M, N, N, Complex.one, A, LDA,
                WORK(IU).asMatrix(LDWRKU), LDWRKU, Complex.zero, U, LDU);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + N;

            // Compute A=Q*R, copying result to U
            // (CWorkspace: need 2*N, prefer N+N*NB)
            // (RWorkspace: 0)

            zgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('L', M, N, A, LDA, U, LDU);

            // Generate Q in U
            // (CWorkspace: need 2*N, prefer N+N*NB)
            // (RWorkspace: 0)

            zungqr(M, N, N, U, LDU, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);

            // Copy R to VT, zeroing out below it

            zlacpy('U', N, N, A, LDA, VT, LDVT);
            if (N > 1) {
              zlaset('L', N - 1, N - 1, Complex.zero, Complex.zero, VT(2, 1),
                  LDVT);
            }
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Bidiagonalize R in VT
            // (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
            // (RWorkspace: need N)

            zgebrd(N, N, VT, LDVT, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply Q in U by left bidiagonalizing vectors
            // in VT
            // (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
            // (RWorkspace: 0)

            zunmbr('Q', 'R', 'N', M, N, N, VT, LDVT, WORK(ITAUQ), U, LDU,
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate right bidiagonalizing vectors in VT
            // (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
            // (RWorkspace: 0)

            zungbr('P', N, N, N, VT, LDVT, WORK(ITAUP), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IRWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of A in U and computing right
            // singular vectors of A in VT
            // (CWorkspace: 0)
            // (RWorkspace: need BDSPAC)

            zbdsqr('U', N, N, M, 0, S, RWORK(IE), VT, LDVT, U, LDU,
                CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
          }
        }
      } else if (WNTUA) {
        if (WNTVN) {
          // Path 7 (M much larger than N, JOBU='A', JOBVT='N')
          // M left singular vectors to be computed in U and
          // no right singular vectors to be computed

          if (LWORK >= N * N + max(N + M, 3 * N)) {
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
            // (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
            // (RWorkspace: 0)

            zgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('L', M, N, A, LDA, U, LDU);

            // Copy R to WORK(IR), zeroing out below it

            zlacpy('U', N, N, A, LDA, WORK(IR).asMatrix(LDWRKR), LDWRKR);
            zlaset('L', N - 1, N - 1, Complex.zero, Complex.zero,
                WORK(IR + 1).asMatrix(LDWRKR), LDWRKR);

            // Generate Q in U
            // (CWorkspace: need N*N+N+M, prefer N*N+N+M*NB)
            // (RWorkspace: 0)

            zungqr(M, M, N, U, LDU, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Bidiagonalize R in WORK(IR)
            // (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
            // (RWorkspace: need N)

            zgebrd(N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, S, RWORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate left bidiagonalizing vectors in WORK(IR)
            // (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
            // (RWorkspace: 0)

            zungbr('Q', N, N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUQ),
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IRWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of R in WORK(IR)
            // (CWorkspace: need N*N)
            // (RWorkspace: need BDSPAC)

            zbdsqr(
                'U',
                N,
                0,
                N,
                0,
                S,
                RWORK(IE),
                CDUM.asMatrix(1),
                1,
                WORK(IR).asMatrix(LDWRKR),
                LDWRKR,
                CDUM.asMatrix(1),
                1,
                RWORK(IRWORK),
                INFO);

            // Multiply Q in U by left singular vectors of R in
            // WORK(IR), storing result in A
            // (CWorkspace: need N*N)
            // (RWorkspace: 0)

            zgemm('N', 'N', M, N, N, Complex.one, U, LDU,
                WORK(IR).asMatrix(LDWRKR), LDWRKR, Complex.zero, A, LDA);

            // Copy left singular vectors of A from A to U

            zlacpy('F', M, N, A, LDA, U, LDU);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + N;

            // Compute A=Q*R, copying result to U
            // (CWorkspace: need 2*N, prefer N+N*NB)
            // (RWorkspace: 0)

            zgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('L', M, N, A, LDA, U, LDU);

            // Generate Q in U
            // (CWorkspace: need N+M, prefer N+M*NB)
            // (RWorkspace: 0)

            zungqr(M, M, N, U, LDU, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Zero out below R in A

            if (N > 1) {
              zlaset(
                  'L', N - 1, N - 1, Complex.zero, Complex.zero, A(2, 1), LDA);
            }

            // Bidiagonalize R in A
            // (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
            // (RWorkspace: need N)

            zgebrd(N, N, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply Q in U by left bidiagonalizing vectors
            // in A
            // (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
            // (RWorkspace: 0)

            zunmbr('Q', 'R', 'N', M, N, N, A, LDA, WORK(ITAUQ), U, LDU,
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IRWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of A in U
            // (CWorkspace: 0)
            // (RWorkspace: need BDSPAC)

            zbdsqr('U', N, 0, M, 0, S, RWORK(IE), CDUM.asMatrix(1), 1, U, LDU,
                CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
          }
        } else if (WNTVO) {
          // Path 8 (M much larger than N, JOBU='A', JOBVT='O')
          // M left singular vectors to be computed in U and
          // N right singular vectors to be overwritten on A

          if (LWORK >= 2 * N * N + max(N + M, 3 * N)) {
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
            // (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)
            // (RWorkspace: 0)

            zgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('L', M, N, A, LDA, U, LDU);

            // Generate Q in U
            // (CWorkspace: need 2*N*N+N+M, prefer 2*N*N+N+M*NB)
            // (RWorkspace: 0)

            zungqr(M, M, N, U, LDU, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);

            // Copy R to WORK(IU), zeroing out below it

            zlacpy('U', N, N, A, LDA, WORK(IU).asMatrix(LDWRKU), LDWRKU);
            zlaset('L', N - 1, N - 1, Complex.zero, Complex.zero,
                WORK(IU + 1).asMatrix(LDWRKU), LDWRKU);
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Bidiagonalize R in WORK(IU), copying result to
            // WORK(IR)
            // (CWorkspace: need   2*N*N+3*N,
            //              prefer 2*N*N+2*N+2*N*NB)
            // (RWorkspace: need   N)

            zgebrd(N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, S, RWORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('U', N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU,
                WORK(IR).asMatrix(LDWRKR), LDWRKR);

            // Generate left bidiagonalizing vectors in WORK(IU)
            // (CWorkspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB)
            // (RWorkspace: 0)

            zungbr('Q', N, N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, WORK(ITAUQ),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate right bidiagonalizing vectors in WORK(IR)
            // (CWorkspace: need   2*N*N+3*N-1,
            //              prefer 2*N*N+2*N+(N-1)*NB)
            // (RWorkspace: 0)

            zungbr('P', N, N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IRWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of R in WORK(IU) and computing
            // right singular vectors of R in WORK(IR)
            // (CWorkspace: need 2*N*N)
            // (RWorkspace: need BDSPAC)

            zbdsqr(
                'U',
                N,
                N,
                N,
                0,
                S,
                RWORK(IE),
                WORK(IR).asMatrix(LDWRKR),
                LDWRKR,
                WORK(IU).asMatrix(LDWRKU),
                LDWRKU,
                CDUM.asMatrix(1),
                1,
                RWORK(IRWORK),
                INFO);

            // Multiply Q in U by left singular vectors of R in
            // WORK(IU), storing result in A
            // (CWorkspace: need N*N)
            // (RWorkspace: 0)

            zgemm('N', 'N', M, N, N, Complex.one, U, LDU,
                WORK(IU).asMatrix(LDWRKU), LDWRKU, Complex.zero, A, LDA);

            // Copy left singular vectors of A from A to U

            zlacpy('F', M, N, A, LDA, U, LDU);

            // Copy right singular vectors of R from WORK(IR) to A

            zlacpy('F', N, N, WORK(IR).asMatrix(LDWRKR), LDWRKR, A, LDA);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + N;

            // Compute A=Q*R, copying result to U
            // (CWorkspace: need 2*N, prefer N+N*NB)
            // (RWorkspace: 0)

            zgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('L', M, N, A, LDA, U, LDU);

            // Generate Q in U
            // (CWorkspace: need N+M, prefer N+M*NB)
            // (RWorkspace: 0)

            zungqr(M, M, N, U, LDU, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Zero out below R in A

            if (N > 1) {
              zlaset(
                  'L', N - 1, N - 1, Complex.zero, Complex.zero, A(2, 1), LDA);
            }

            // Bidiagonalize R in A
            // (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
            // (RWorkspace: need N)

            zgebrd(N, N, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply Q in U by left bidiagonalizing vectors
            // in A
            // (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
            // (RWorkspace: 0)

            zunmbr('Q', 'R', 'N', M, N, N, A, LDA, WORK(ITAUQ), U, LDU,
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate right bidiagonalizing vectors in A
            // (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
            // (RWorkspace: 0)

            zungbr('P', N, N, N, A, LDA, WORK(ITAUP), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IRWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of A in U and computing right
            // singular vectors of A in A
            // (CWorkspace: 0)
            // (RWorkspace: need BDSPAC)

            zbdsqr('U', N, N, M, 0, S, RWORK(IE), A, LDA, U, LDU,
                CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
          }
        } else if (WNTVAS) {
          // Path 9 (M much larger than N, JOBU='A', JOBVT='S'
          //         or 'A')
          // M left singular vectors to be computed in U and
          // N right singular vectors to be computed in VT

          if (LWORK >= N * N + max(N + M, 3 * N)) {
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
            // (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
            // (RWorkspace: 0)

            zgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('L', M, N, A, LDA, U, LDU);

            // Generate Q in U
            // (CWorkspace: need N*N+N+M, prefer N*N+N+M*NB)
            // (RWorkspace: 0)

            zungqr(M, M, N, U, LDU, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);

            // Copy R to WORK(IU), zeroing out below it

            zlacpy('U', N, N, A, LDA, WORK(IU).asMatrix(LDWRKU), LDWRKU);
            zlaset('L', N - 1, N - 1, Complex.zero, Complex.zero,
                WORK(IU + 1).asMatrix(LDWRKU), LDWRKU);
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Bidiagonalize R in WORK(IU), copying result to VT
            // (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
            // (RWorkspace: need N)

            zgebrd(N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, S, RWORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('U', N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, VT, LDVT);

            // Generate left bidiagonalizing vectors in WORK(IU)
            // (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
            // (RWorkspace: 0)

            zungbr('Q', N, N, N, WORK(IU).asMatrix(LDWRKU), LDWRKU, WORK(ITAUQ),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate right bidiagonalizing vectors in VT
            // (CWorkspace: need   N*N+3*N-1,
            //              prefer N*N+2*N+(N-1)*NB)
            // (RWorkspace: need   0)

            zungbr('P', N, N, N, VT, LDVT, WORK(ITAUP), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IRWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of R in WORK(IU) and computing
            // right singular vectors of R in VT
            // (CWorkspace: need N*N)
            // (RWorkspace: need BDSPAC)

            zbdsqr(
                'U',
                N,
                N,
                N,
                0,
                S,
                RWORK(IE),
                VT,
                LDVT,
                WORK(IU).asMatrix(LDWRKU),
                LDWRKU,
                CDUM.asMatrix(1),
                1,
                RWORK(IRWORK),
                INFO);

            // Multiply Q in U by left singular vectors of R in
            // WORK(IU), storing result in A
            // (CWorkspace: need N*N)
            // (RWorkspace: 0)

            zgemm('N', 'N', M, N, N, Complex.one, U, LDU,
                WORK(IU).asMatrix(LDWRKU), LDWRKU, Complex.zero, A, LDA);

            // Copy left singular vectors of A from A to U

            zlacpy('F', M, N, A, LDA, U, LDU);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + N;

            // Compute A=Q*R, copying result to U
            // (CWorkspace: need 2*N, prefer N+N*NB)
            // (RWorkspace: 0)

            zgeqrf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('L', M, N, A, LDA, U, LDU);

            // Generate Q in U
            // (CWorkspace: need N+M, prefer N+M*NB)
            // (RWorkspace: 0)

            zungqr(M, M, N, U, LDU, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);

            // Copy R from A to VT, zeroing out below it

            zlacpy('U', N, N, A, LDA, VT, LDVT);
            if (N > 1) {
              zlaset('L', N - 1, N - 1, Complex.zero, Complex.zero, VT(2, 1),
                  LDVT);
            }
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Bidiagonalize R in VT
            // (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
            // (RWorkspace: need N)

            zgebrd(N, N, VT, LDVT, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply Q in U by left bidiagonalizing vectors
            // in VT
            // (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
            // (RWorkspace: 0)

            zunmbr('Q', 'R', 'N', M, N, N, VT, LDVT, WORK(ITAUQ), U, LDU,
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate right bidiagonalizing vectors in VT
            // (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
            // (RWorkspace: 0)

            zungbr('P', N, N, N, VT, LDVT, WORK(ITAUP), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IRWORK = IE + N;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of A in U and computing right
            // singular vectors of A in VT
            // (CWorkspace: 0)
            // (RWorkspace: need BDSPAC)

            zbdsqr('U', N, N, M, 0, S, RWORK(IE), VT, LDVT, U, LDU,
                CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
          }
        }
      }
    } else {
      // M < MNTHR

      // Path 10 (M at least N, but not much larger)
      // Reduce to bidiagonal form without QR decomposition

      IE = 1;
      ITAUQ = 1;
      ITAUP = ITAUQ + N;
      IWORK = ITAUP + N;

      // Bidiagonalize A
      // (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB)
      // (RWorkspace: need N)

      zgebrd(M, N, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP), WORK(IWORK),
          LWORK - IWORK + 1, IERR);
      if (WNTUAS) {
        // If left singular vectors desired in U, copy result to U
        // and generate left bidiagonalizing vectors in U
        // (CWorkspace: need 2*N+NCU, prefer 2*N+NCU*NB)
        // (RWorkspace: 0)

        zlacpy('L', M, N, A, LDA, U, LDU);
        if (WNTUS) NCU = N;
        if (WNTUA) NCU = M;
        zungbr('Q', M, NCU, N, U, LDU, WORK(ITAUQ), WORK(IWORK),
            LWORK - IWORK + 1, IERR);
      }
      if (WNTVAS) {
        // If right singular vectors desired in VT, copy result to
        // VT and generate right bidiagonalizing vectors in VT
        // (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
        // (RWorkspace: 0)

        zlacpy('U', N, N, A, LDA, VT, LDVT);
        zungbr('P', N, N, N, VT, LDVT, WORK(ITAUP), WORK(IWORK),
            LWORK - IWORK + 1, IERR);
      }
      if (WNTUO) {
        // If left singular vectors desired in A, generate left
        // bidiagonalizing vectors in A
        // (CWorkspace: need 3*N, prefer 2*N+N*NB)
        // (RWorkspace: 0)

        zungbr('Q', M, N, N, A, LDA, WORK(ITAUQ), WORK(IWORK),
            LWORK - IWORK + 1, IERR);
      }
      if (WNTVO) {
        // If right singular vectors desired in A, generate right
        // bidiagonalizing vectors in A
        // (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
        // (RWorkspace: 0)

        zungbr('P', N, N, N, A, LDA, WORK(ITAUP), WORK(IWORK),
            LWORK - IWORK + 1, IERR);
      }
      IRWORK = IE + N;
      if (WNTUAS || WNTUO) NRU = M;
      if (WNTUN) NRU = 0;
      if (WNTVAS || WNTVO) NCVT = N;
      if (WNTVN) NCVT = 0;
      if (!WNTUO && !WNTVO) {
        // Perform bidiagonal QR iteration, if desired, computing
        // left singular vectors in U and computing right singular
        // vectors in VT
        // (CWorkspace: 0)
        // (RWorkspace: need BDSPAC)

        zbdsqr('U', N, NCVT, NRU, 0, S, RWORK(IE), VT, LDVT, U, LDU,
            CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
      } else if (!WNTUO && WNTVO) {
        // Perform bidiagonal QR iteration, if desired, computing
        // left singular vectors in U and computing right singular
        // vectors in A
        // (CWorkspace: 0)
        // (RWorkspace: need BDSPAC)

        zbdsqr('U', N, NCVT, NRU, 0, S, RWORK(IE), A, LDA, U, LDU,
            CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
      } else {
        // Perform bidiagonal QR iteration, if desired, computing
        // left singular vectors in A and computing right singular
        // vectors in VT
        // (CWorkspace: 0)
        // (RWorkspace: need BDSPAC)

        zbdsqr('U', N, NCVT, NRU, 0, S, RWORK(IE), VT, LDVT, A, LDA,
            CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
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
        // (CWorkspace: need 2*M, prefer M+M*NB)
        // (RWorkspace: 0)

        zgelqf(M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

        // Zero out above L

        zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero, A(1, 2), LDA);
        IE = 1;
        ITAUQ = 1;
        ITAUP = ITAUQ + M;
        IWORK = ITAUP + M;

        // Bidiagonalize L in A
        // (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
        // (RWorkspace: need M)

        zgebrd(M, M, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
            WORK(IWORK), LWORK - IWORK + 1, IERR);
        if (WNTUO || WNTUAS) {
          // If left singular vectors desired, generate Q
          // (CWorkspace: need 3*M, prefer 2*M+M*NB)
          // (RWorkspace: 0)

          zungbr('Q', M, M, M, A, LDA, WORK(ITAUQ), WORK(IWORK),
              LWORK - IWORK + 1, IERR);
        }
        IRWORK = IE + M;
        NRU = 0;
        if (WNTUO || WNTUAS) NRU = M;

        // Perform bidiagonal QR iteration, computing left singular
        // vectors of A in A if desired
        // (CWorkspace: 0)
        // (RWorkspace: need BDSPAC)

        zbdsqr('U', M, 0, NRU, 0, S, RWORK(IE), CDUM.asMatrix(1), 1, A, LDA,
            CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);

        // If left singular vectors desired in U, copy them there

        if (WNTUAS) zlacpy('F', M, M, A, LDA, U, LDU);
      } else if (WNTVO && WNTUN) {
        // Path 2t(N much larger than M, JOBU='N', JOBVT='O')
        // M right singular vectors to be overwritten on A and
        // no left singular vectors to be computed

        if (LWORK >= M * M + 3 * M) {
          // Sufficient workspace for a fast algorithm

          IR = 1;
          if (LWORK >= max(WRKBL, LDA * N) + LDA * M) {
            // WORK(IU) is LDA by N and WORK(IR) is LDA by M

            LDWRKU = LDA;
            CHUNK = N;
            LDWRKR = LDA;
          } else if (LWORK >= max(WRKBL, LDA * N) + M * M) {
            // WORK(IU) is LDA by N and WORK(IR) is M by M

            LDWRKU = LDA;
            CHUNK = N;
            LDWRKR = M;
          } else {
            // WORK(IU) is M by CHUNK and WORK(IR) is M by M

            LDWRKU = M;
            CHUNK = (LWORK - M * M) ~/ M;
            LDWRKR = M;
          }
          ITAU = IR + LDWRKR * M;
          IWORK = ITAU + M;

          // Compute A=L*Q
          // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
          // (RWorkspace: 0)

          zgelqf(
              M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Copy L to WORK(IR) and zero out above it

          zlacpy('L', M, M, A, LDA, WORK(IR).asMatrix(LDWRKR), LDWRKR);
          zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero,
              WORK(IR + LDWRKR).asMatrix(LDWRKR), LDWRKR);

          // Generate Q in A
          // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
          // (RWorkspace: 0)

          zunglq(M, N, M, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
              IERR);
          IE = 1;
          ITAUQ = ITAU;
          ITAUP = ITAUQ + M;
          IWORK = ITAUP + M;

          // Bidiagonalize L in WORK(IR)
          // (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
          // (RWorkspace: need M)

          zgebrd(M, M, WORK(IR).asMatrix(LDWRKR), LDWRKR, S, RWORK(IE),
              WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Generate right vectors bidiagonalizing L
          // (CWorkspace: need M*M+3*M-1, prefer M*M+2*M+(M-1)*NB)
          // (RWorkspace: 0)

          zungbr('P', M, M, M, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUP),
              WORK(IWORK), LWORK - IWORK + 1, IERR);
          IRWORK = IE + M;

          // Perform bidiagonal QR iteration, computing right
          // singular vectors of L in WORK(IR)
          // (CWorkspace: need M*M)
          // (RWorkspace: need BDSPAC)

          zbdsqr(
              'U',
              M,
              M,
              0,
              0,
              S,
              RWORK(IE),
              WORK(IR).asMatrix(LDWRKR),
              LDWRKR,
              CDUM.asMatrix(1),
              1,
              CDUM.asMatrix(1),
              1,
              RWORK(IRWORK),
              INFO);
          IU = ITAUQ;

          // Multiply right singular vectors of L in WORK(IR) by Q
          // in A, storing result in WORK(IU) and copying to A
          // (CWorkspace: need M*M+M, prefer M*M+M*N)
          // (RWorkspace: 0)

          for (I = 1; I <= N; I += CHUNK) {
            BLK = min(N - I + 1, CHUNK);
            zgemm(
                'N',
                'N',
                M,
                BLK,
                M,
                Complex.one,
                WORK(IR).asMatrix(LDWRKR),
                LDWRKR,
                A(1, I),
                LDA,
                Complex.zero,
                WORK(IU).asMatrix(LDWRKU),
                LDWRKU);
            zlacpy(
                'F', M, BLK, WORK(IU).asMatrix(LDWRKU), LDWRKU, A(1, I), LDA);
          }
        } else {
          // Insufficient workspace for a fast algorithm

          IE = 1;
          ITAUQ = 1;
          ITAUP = ITAUQ + M;
          IWORK = ITAUP + M;

          // Bidiagonalize A
          // (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB)
          // (RWorkspace: need M)

          zgebrd(M, N, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
              WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Generate right vectors bidiagonalizing A
          // (CWorkspace: need 3*M, prefer 2*M+M*NB)
          // (RWorkspace: 0)

          zungbr('P', M, N, M, A, LDA, WORK(ITAUP), WORK(IWORK),
              LWORK - IWORK + 1, IERR);
          IRWORK = IE + M;

          // Perform bidiagonal QR iteration, computing right
          // singular vectors of A in A
          // (CWorkspace: 0)
          // (RWorkspace: need BDSPAC)

          zbdsqr('L', M, N, 0, 0, S, RWORK(IE), A, LDA, CDUM.asMatrix(1), 1,
              CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
        }
      } else if (WNTVO && WNTUAS) {
        // Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O')
        // M right singular vectors to be overwritten on A and
        // M left singular vectors to be computed in U

        if (LWORK >= M * M + 3 * M) {
          // Sufficient workspace for a fast algorithm

          IR = 1;
          if (LWORK >= max(WRKBL, LDA * N) + LDA * M) {
            // WORK(IU) is LDA by N and WORK(IR) is LDA by M

            LDWRKU = LDA;
            CHUNK = N;
            LDWRKR = LDA;
          } else if (LWORK >= max(WRKBL, LDA * N) + M * M) {
            // WORK(IU) is LDA by N and WORK(IR) is M by M

            LDWRKU = LDA;
            CHUNK = N;
            LDWRKR = M;
          } else {
            // WORK(IU) is M by CHUNK and WORK(IR) is M by M

            LDWRKU = M;
            CHUNK = (LWORK - M * M) ~/ M;
            LDWRKR = M;
          }
          ITAU = IR + LDWRKR * M;
          IWORK = ITAU + M;

          // Compute A=L*Q
          // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
          // (RWorkspace: 0)

          zgelqf(
              M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Copy L to U, zeroing about above it

          zlacpy('L', M, M, A, LDA, U, LDU);
          zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero, U(1, 2), LDU);

          // Generate Q in A
          // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
          // (RWorkspace: 0)

          zunglq(M, N, M, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
              IERR);
          IE = 1;
          ITAUQ = ITAU;
          ITAUP = ITAUQ + M;
          IWORK = ITAUP + M;

          // Bidiagonalize L in U, copying result to WORK(IR)
          // (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
          // (RWorkspace: need M)

          zgebrd(M, M, U, LDU, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
              WORK(IWORK), LWORK - IWORK + 1, IERR);
          zlacpy('U', M, M, U, LDU, WORK(IR).asMatrix(LDWRKR), LDWRKR);

          // Generate right vectors bidiagonalizing L in WORK(IR)
          // (CWorkspace: need M*M+3*M-1, prefer M*M+2*M+(M-1)*NB)
          // (RWorkspace: 0)

          zungbr('P', M, M, M, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUP),
              WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Generate left vectors bidiagonalizing L in U
          // (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB)
          // (RWorkspace: 0)

          zungbr('Q', M, M, M, U, LDU, WORK(ITAUQ), WORK(IWORK),
              LWORK - IWORK + 1, IERR);
          IRWORK = IE + M;

          // Perform bidiagonal QR iteration, computing left
          // singular vectors of L in U, and computing right
          // singular vectors of L in WORK(IR)
          // (CWorkspace: need M*M)
          // (RWorkspace: need BDSPAC)

          zbdsqr('U', M, M, M, 0, S, RWORK(IE), WORK(IR).asMatrix(LDWRKR),
              LDWRKR, U, LDU, CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
          IU = ITAUQ;

          // Multiply right singular vectors of L in WORK(IR) by Q
          // in A, storing result in WORK(IU) and copying to A
          // (CWorkspace: need M*M+M, prefer M*M+M*N))
          // (RWorkspace: 0)

          for (I = 1; I <= N; I += CHUNK) {
            BLK = min(N - I + 1, CHUNK);
            zgemm(
                'N',
                'N',
                M,
                BLK,
                M,
                Complex.one,
                WORK(IR).asMatrix(LDWRKR),
                LDWRKR,
                A(1, I),
                LDA,
                Complex.zero,
                WORK(IU).asMatrix(LDWRKU),
                LDWRKU);
            zlacpy(
                'F', M, BLK, WORK(IU).asMatrix(LDWRKU), LDWRKU, A(1, I), LDA);
          }
        } else {
          // Insufficient workspace for a fast algorithm

          ITAU = 1;
          IWORK = ITAU + M;

          // Compute A=L*Q
          // (CWorkspace: need 2*M, prefer M+M*NB)
          // (RWorkspace: 0)

          zgelqf(
              M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Copy L to U, zeroing out above it

          zlacpy('L', M, M, A, LDA, U, LDU);
          zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero, U(1, 2), LDU);

          // Generate Q in A
          // (CWorkspace: need 2*M, prefer M+M*NB)
          // (RWorkspace: 0)

          zunglq(M, N, M, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
              IERR);
          IE = 1;
          ITAUQ = ITAU;
          ITAUP = ITAUQ + M;
          IWORK = ITAUP + M;

          // Bidiagonalize L in U
          // (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
          // (RWorkspace: need M)

          zgebrd(M, M, U, LDU, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
              WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Multiply right vectors bidiagonalizing L by Q in A
          // (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
          // (RWorkspace: 0)

          zunmbr('P', 'L', 'C', M, N, M, U, LDU, WORK(ITAUP), A, LDA,
              WORK(IWORK), LWORK - IWORK + 1, IERR);

          // Generate left vectors bidiagonalizing L in U
          // (CWorkspace: need 3*M, prefer 2*M+M*NB)
          // (RWorkspace: 0)

          zungbr('Q', M, M, M, U, LDU, WORK(ITAUQ), WORK(IWORK),
              LWORK - IWORK + 1, IERR);
          IRWORK = IE + M;

          // Perform bidiagonal QR iteration, computing left
          // singular vectors of A in U and computing right
          // singular vectors of A in A
          // (CWorkspace: 0)
          // (RWorkspace: need BDSPAC)

          zbdsqr('U', M, N, M, 0, S, RWORK(IE), A, LDA, U, LDU,
              CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
        }
      } else if (WNTVS) {
        if (WNTUN) {
          // Path 4t(N much larger than M, JOBU='N', JOBVT='S')
          // M right singular vectors to be computed in VT and
          // no left singular vectors to be computed

          if (LWORK >= M * M + 3 * M) {
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
            // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
            // (RWorkspace: 0)

            zgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Copy L to WORK(IR), zeroing out above it

            zlacpy('L', M, M, A, LDA, WORK(IR).asMatrix(LDWRKR), LDWRKR);
            zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero,
                WORK(IR + LDWRKR).asMatrix(LDWRKR), LDWRKR);

            // Generate Q in A
            // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
            // (RWorkspace: 0)

            zunglq(M, N, M, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Bidiagonalize L in WORK(IR)
            // (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
            // (RWorkspace: need M)

            zgebrd(M, M, WORK(IR).asMatrix(LDWRKR), LDWRKR, S, RWORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate right vectors bidiagonalizing L in
            // WORK(IR)
            // (CWorkspace: need M*M+3*M, prefer M*M+2*M+(M-1)*NB)
            // (RWorkspace: 0)

            zungbr('P', M, M, M, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IRWORK = IE + M;

            // Perform bidiagonal QR iteration, computing right
            // singular vectors of L in WORK(IR)
            // (CWorkspace: need M*M)
            // (RWorkspace: need BDSPAC)

            zbdsqr(
                'U',
                M,
                M,
                0,
                0,
                S,
                RWORK(IE),
                WORK(IR).asMatrix(LDWRKR),
                LDWRKR,
                CDUM.asMatrix(1),
                1,
                CDUM.asMatrix(1),
                1,
                RWORK(IRWORK),
                INFO);

            // Multiply right singular vectors of L in WORK(IR) by
            // Q in A, storing result in VT
            // (CWorkspace: need M*M)
            // (RWorkspace: 0)

            zgemm('N', 'N', M, N, M, Complex.one, WORK(IR).asMatrix(LDWRKR),
                LDWRKR, A, LDA, Complex.zero, VT, LDVT);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + M;

            // Compute A=L*Q
            // (CWorkspace: need 2*M, prefer M+M*NB)
            // (RWorkspace: 0)

            zgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Copy result to VT

            zlacpy('U', M, N, A, LDA, VT, LDVT);

            // Generate Q in VT
            // (CWorkspace: need 2*M, prefer M+M*NB)
            // (RWorkspace: 0)

            zunglq(M, N, M, VT, LDVT, WORK(ITAU), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Zero out above L in A

            zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero, A(1, 2), LDA);

            // Bidiagonalize L in A
            // (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
            // (RWorkspace: need M)

            zgebrd(M, M, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply right vectors bidiagonalizing L by Q in VT
            // (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
            // (RWorkspace: 0)

            zunmbr('P', 'L', 'C', M, N, M, A, LDA, WORK(ITAUP), VT, LDVT,
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IRWORK = IE + M;

            // Perform bidiagonal QR iteration, computing right
            // singular vectors of A in VT
            // (CWorkspace: 0)
            // (RWorkspace: need BDSPAC)

            zbdsqr('U', M, N, 0, 0, S, RWORK(IE), VT, LDVT, CDUM.asMatrix(1), 1,
                CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
          }
        } else if (WNTUO) {
          // Path 5t(N much larger than M, JOBU='O', JOBVT='S')
          // M right singular vectors to be computed in VT and
          // M left singular vectors to be overwritten on A

          if (LWORK >= 2 * M * M + 3 * M) {
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
            // (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)
            // (RWorkspace: 0)

            zgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Copy L to WORK(IU), zeroing out below it

            zlacpy('L', M, M, A, LDA, WORK(IU).asMatrix(LDWRKU), LDWRKU);
            zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero,
                WORK(IU + LDWRKU).asMatrix(LDWRKU), LDWRKU);

            // Generate Q in A
            // (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)
            // (RWorkspace: 0)

            zunglq(M, N, M, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Bidiagonalize L in WORK(IU), copying result to
            // WORK(IR)
            // (CWorkspace: need   2*M*M+3*M,
            //              prefer 2*M*M+2*M+2*M*NB)
            // (RWorkspace: need   M)

            zgebrd(M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU, S, RWORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('L', M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU,
                WORK(IR).asMatrix(LDWRKR), LDWRKR);

            // Generate right bidiagonalizing vectors in WORK(IU)
            // (CWorkspace: need   2*M*M+3*M-1,
            //              prefer 2*M*M+2*M+(M-1)*NB)
            // (RWorkspace: 0)

            zungbr('P', M, M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU, WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate left bidiagonalizing vectors in WORK(IR)
            // (CWorkspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB)
            // (RWorkspace: 0)

            zungbr('Q', M, M, M, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUQ),
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IRWORK = IE + M;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of L in WORK(IR) and computing
            // right singular vectors of L in WORK(IU)
            // (CWorkspace: need 2*M*M)
            // (RWorkspace: need BDSPAC)

            zbdsqr(
                'U',
                M,
                M,
                M,
                0,
                S,
                RWORK(IE),
                WORK(IU).asMatrix(LDWRKU),
                LDWRKU,
                WORK(IR).asMatrix(LDWRKR),
                LDWRKR,
                CDUM.asMatrix(1),
                1,
                RWORK(IRWORK),
                INFO);

            // Multiply right singular vectors of L in WORK(IU) by
            // Q in A, storing result in VT
            // (CWorkspace: need M*M)
            // (RWorkspace: 0)

            zgemm('N', 'N', M, N, M, Complex.one, WORK(IU).asMatrix(LDWRKU),
                LDWRKU, A, LDA, Complex.zero, VT, LDVT);

            // Copy left singular vectors of L to A
            // (CWorkspace: need M*M)
            // (RWorkspace: 0)

            zlacpy('F', M, M, WORK(IR).asMatrix(LDWRKR), LDWRKR, A, LDA);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + M;

            // Compute A=L*Q, copying result to VT
            // (CWorkspace: need 2*M, prefer M+M*NB)
            // (RWorkspace: 0)

            zgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('U', M, N, A, LDA, VT, LDVT);

            // Generate Q in VT
            // (CWorkspace: need 2*M, prefer M+M*NB)
            // (RWorkspace: 0)

            zunglq(M, N, M, VT, LDVT, WORK(ITAU), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Zero out above L in A

            zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero, A(1, 2), LDA);

            // Bidiagonalize L in A
            // (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
            // (RWorkspace: need M)

            zgebrd(M, M, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply right vectors bidiagonalizing L by Q in VT
            // (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
            // (RWorkspace: 0)

            zunmbr('P', 'L', 'C', M, N, M, A, LDA, WORK(ITAUP), VT, LDVT,
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate left bidiagonalizing vectors of L in A
            // (CWorkspace: need 3*M, prefer 2*M+M*NB)
            // (RWorkspace: 0)

            zungbr('Q', M, M, M, A, LDA, WORK(ITAUQ), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IRWORK = IE + M;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of A in A and computing right
            // singular vectors of A in VT
            // (CWorkspace: 0)
            // (RWorkspace: need BDSPAC)

            zbdsqr('U', M, N, M, 0, S, RWORK(IE), VT, LDVT, A, LDA,
                CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
          }
        } else if (WNTUAS) {
          // Path 6t(N much larger than M, JOBU='S' or 'A',
          //         JOBVT='S')
          // M right singular vectors to be computed in VT and
          // M left singular vectors to be computed in U

          if (LWORK >= M * M + 3 * M) {
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
            // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
            // (RWorkspace: 0)

            zgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Copy L to WORK(IU), zeroing out above it

            zlacpy('L', M, M, A, LDA, WORK(IU).asMatrix(LDWRKU), LDWRKU);
            zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero,
                WORK(IU + LDWRKU).asMatrix(LDWRKU), LDWRKU);

            // Generate Q in A
            // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
            // (RWorkspace: 0)

            zunglq(M, N, M, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1,
                IERR);
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Bidiagonalize L in WORK(IU), copying result to U
            // (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
            // (RWorkspace: need M)

            zgebrd(M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU, S, RWORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('L', M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU, U, LDU);

            // Generate right bidiagonalizing vectors in WORK(IU)
            // (CWorkspace: need   M*M+3*M-1,
            //              prefer M*M+2*M+(M-1)*NB)
            // (RWorkspace: 0)

            zungbr('P', M, M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU, WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate left bidiagonalizing vectors in U
            // (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB)
            // (RWorkspace: 0)

            zungbr('Q', M, M, M, U, LDU, WORK(ITAUQ), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IRWORK = IE + M;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of L in U and computing right
            // singular vectors of L in WORK(IU)
            // (CWorkspace: need M*M)
            // (RWorkspace: need BDSPAC)

            zbdsqr('U', M, M, M, 0, S, RWORK(IE), WORK(IU).asMatrix(LDWRKU),
                LDWRKU, U, LDU, CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);

            // Multiply right singular vectors of L in WORK(IU) by
            // Q in A, storing result in VT
            // (CWorkspace: need M*M)
            // (RWorkspace: 0)

            zgemm('N', 'N', M, N, M, Complex.one, WORK(IU).asMatrix(LDWRKU),
                LDWRKU, A, LDA, Complex.zero, VT, LDVT);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + M;

            // Compute A=L*Q, copying result to VT
            // (CWorkspace: need 2*M, prefer M+M*NB)
            // (RWorkspace: 0)

            zgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('U', M, N, A, LDA, VT, LDVT);

            // Generate Q in VT
            // (CWorkspace: need 2*M, prefer M+M*NB)
            // (RWorkspace: 0)

            zunglq(M, N, M, VT, LDVT, WORK(ITAU), WORK(IWORK),
                LWORK - IWORK + 1, IERR);

            // Copy L to U, zeroing out above it

            zlacpy('L', M, M, A, LDA, U, LDU);
            zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero, U(1, 2), LDU);
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Bidiagonalize L in U
            // (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
            // (RWorkspace: need M)

            zgebrd(M, M, U, LDU, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply right bidiagonalizing vectors in U by Q
            // in VT
            // (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
            // (RWorkspace: 0)

            zunmbr('P', 'L', 'C', M, N, M, U, LDU, WORK(ITAUP), VT, LDVT,
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate left bidiagonalizing vectors in U
            // (CWorkspace: need 3*M, prefer 2*M+M*NB)
            // (RWorkspace: 0)

            zungbr('Q', M, M, M, U, LDU, WORK(ITAUQ), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IRWORK = IE + M;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of A in U and computing right
            // singular vectors of A in VT
            // (CWorkspace: 0)
            // (RWorkspace: need BDSPAC)

            zbdsqr('U', M, N, M, 0, S, RWORK(IE), VT, LDVT, U, LDU,
                CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
          }
        }
      } else if (WNTVA) {
        if (WNTUN) {
          // Path 7t(N much larger than M, JOBU='N', JOBVT='A')
          // N right singular vectors to be computed in VT and
          // no left singular vectors to be computed

          if (LWORK >= M * M + max(N + M, 3 * M)) {
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
            // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
            // (RWorkspace: 0)

            zgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('U', M, N, A, LDA, VT, LDVT);

            // Copy L to WORK(IR), zeroing out above it

            zlacpy('L', M, M, A, LDA, WORK(IR).asMatrix(LDWRKR), LDWRKR);
            zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero,
                WORK(IR + LDWRKR).asMatrix(LDWRKR), LDWRKR);

            // Generate Q in VT
            // (CWorkspace: need M*M+M+N, prefer M*M+M+N*NB)
            // (RWorkspace: 0)

            zunglq(N, N, M, VT, LDVT, WORK(ITAU), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Bidiagonalize L in WORK(IR)
            // (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
            // (RWorkspace: need M)

            zgebrd(M, M, WORK(IR).asMatrix(LDWRKR), LDWRKR, S, RWORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate right bidiagonalizing vectors in WORK(IR)
            // (CWorkspace: need   M*M+3*M-1,
            //              prefer M*M+2*M+(M-1)*NB)
            // (RWorkspace: 0)

            zungbr('P', M, M, M, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IRWORK = IE + M;

            // Perform bidiagonal QR iteration, computing right
            // singular vectors of L in WORK(IR)
            // (CWorkspace: need M*M)
            // (RWorkspace: need BDSPAC)

            zbdsqr(
                'U',
                M,
                M,
                0,
                0,
                S,
                RWORK(IE),
                WORK(IR).asMatrix(LDWRKR),
                LDWRKR,
                CDUM.asMatrix(1),
                1,
                CDUM.asMatrix(1),
                1,
                RWORK(IRWORK),
                INFO);

            // Multiply right singular vectors of L in WORK(IR) by
            // Q in VT, storing result in A
            // (CWorkspace: need M*M)
            // (RWorkspace: 0)

            zgemm('N', 'N', M, N, M, Complex.one, WORK(IR).asMatrix(LDWRKR),
                LDWRKR, VT, LDVT, Complex.zero, A, LDA);

            // Copy right singular vectors of A from A to VT

            zlacpy('F', M, N, A, LDA, VT, LDVT);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + M;

            // Compute A=L*Q, copying result to VT
            // (CWorkspace: need 2*M, prefer M+M*NB)
            // (RWorkspace: 0)

            zgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('U', M, N, A, LDA, VT, LDVT);

            // Generate Q in VT
            // (CWorkspace: need M+N, prefer M+N*NB)
            // (RWorkspace: 0)

            zunglq(N, N, M, VT, LDVT, WORK(ITAU), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Zero out above L in A

            zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero, A(1, 2), LDA);

            // Bidiagonalize L in A
            // (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
            // (RWorkspace: need M)

            zgebrd(M, M, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply right bidiagonalizing vectors in A by Q
            // in VT
            // (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
            // (RWorkspace: 0)

            zunmbr('P', 'L', 'C', M, N, M, A, LDA, WORK(ITAUP), VT, LDVT,
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IRWORK = IE + M;

            // Perform bidiagonal QR iteration, computing right
            // singular vectors of A in VT
            // (CWorkspace: 0)
            // (RWorkspace: need BDSPAC)

            zbdsqr('U', M, N, 0, 0, S, RWORK(IE), VT, LDVT, CDUM.asMatrix(1), 1,
                CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
          }
        } else if (WNTUO) {
          // Path 8t(N much larger than M, JOBU='O', JOBVT='A')
          // N right singular vectors to be computed in VT and
          // M left singular vectors to be overwritten on A

          if (LWORK >= 2 * M * M + max(N + M, 3 * M)) {
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
            // (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)
            // (RWorkspace: 0)

            zgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('U', M, N, A, LDA, VT, LDVT);

            // Generate Q in VT
            // (CWorkspace: need 2*M*M+M+N, prefer 2*M*M+M+N*NB)
            // (RWorkspace: 0)

            zunglq(N, N, M, VT, LDVT, WORK(ITAU), WORK(IWORK),
                LWORK - IWORK + 1, IERR);

            // Copy L to WORK(IU), zeroing out above it

            zlacpy('L', M, M, A, LDA, WORK(IU).asMatrix(LDWRKU), LDWRKU);
            zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero,
                WORK(IU + LDWRKU).asMatrix(LDWRKU), LDWRKU);
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Bidiagonalize L in WORK(IU), copying result to
            // WORK(IR)
            // (CWorkspace: need   2*M*M+3*M,
            //              prefer 2*M*M+2*M+2*M*NB)
            // (RWorkspace: need   M)

            zgebrd(M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU, S, RWORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('L', M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU,
                WORK(IR).asMatrix(LDWRKR), LDWRKR);

            // Generate right bidiagonalizing vectors in WORK(IU)
            // (CWorkspace: need   2*M*M+3*M-1,
            //              prefer 2*M*M+2*M+(M-1)*NB)
            // (RWorkspace: 0)

            zungbr('P', M, M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU, WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate left bidiagonalizing vectors in WORK(IR)
            // (CWorkspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB)
            // (RWorkspace: 0)

            zungbr('Q', M, M, M, WORK(IR).asMatrix(LDWRKR), LDWRKR, WORK(ITAUQ),
                WORK(IWORK), LWORK - IWORK + 1, IERR);
            IRWORK = IE + M;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of L in WORK(IR) and computing
            // right singular vectors of L in WORK(IU)
            // (CWorkspace: need 2*M*M)
            // (RWorkspace: need BDSPAC)

            zbdsqr(
                'U',
                M,
                M,
                M,
                0,
                S,
                RWORK(IE),
                WORK(IU).asMatrix(LDWRKU),
                LDWRKU,
                WORK(IR).asMatrix(LDWRKR),
                LDWRKR,
                CDUM.asMatrix(1),
                1,
                RWORK(IRWORK),
                INFO);

            // Multiply right singular vectors of L in WORK(IU) by
            // Q in VT, storing result in A
            // (CWorkspace: need M*M)
            // (RWorkspace: 0)

            zgemm('N', 'N', M, N, M, Complex.one, WORK(IU).asMatrix(LDWRKU),
                LDWRKU, VT, LDVT, Complex.zero, A, LDA);

            // Copy right singular vectors of A from A to VT

            zlacpy('F', M, N, A, LDA, VT, LDVT);

            // Copy left singular vectors of A from WORK(IR) to A

            zlacpy('F', M, M, WORK(IR).asMatrix(LDWRKR), LDWRKR, A, LDA);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + M;

            // Compute A=L*Q, copying result to VT
            // (CWorkspace: need 2*M, prefer M+M*NB)
            // (RWorkspace: 0)

            zgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('U', M, N, A, LDA, VT, LDVT);

            // Generate Q in VT
            // (CWorkspace: need M+N, prefer M+N*NB)
            // (RWorkspace: 0)

            zunglq(N, N, M, VT, LDVT, WORK(ITAU), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Zero out above L in A

            zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero, A(1, 2), LDA);

            // Bidiagonalize L in A
            // (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
            // (RWorkspace: need M)

            zgebrd(M, M, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply right bidiagonalizing vectors in A by Q
            // in VT
            // (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
            // (RWorkspace: 0)

            zunmbr('P', 'L', 'C', M, N, M, A, LDA, WORK(ITAUP), VT, LDVT,
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate left bidiagonalizing vectors in A
            // (CWorkspace: need 3*M, prefer 2*M+M*NB)
            // (RWorkspace: 0)

            zungbr('Q', M, M, M, A, LDA, WORK(ITAUQ), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IRWORK = IE + M;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of A in A and computing right
            // singular vectors of A in VT
            // (CWorkspace: 0)
            // (RWorkspace: need BDSPAC)

            zbdsqr('U', M, N, M, 0, S, RWORK(IE), VT, LDVT, A, LDA,
                CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
          }
        } else if (WNTUAS) {
          // Path 9t(N much larger than M, JOBU='S' or 'A',
          //         JOBVT='A')
          // N right singular vectors to be computed in VT and
          // M left singular vectors to be computed in U

          if (LWORK >= M * M + max(N + M, 3 * M)) {
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
            // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
            // (RWorkspace: 0)

            zgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('U', M, N, A, LDA, VT, LDVT);

            // Generate Q in VT
            // (CWorkspace: need M*M+M+N, prefer M*M+M+N*NB)
            // (RWorkspace: 0)

            zunglq(N, N, M, VT, LDVT, WORK(ITAU), WORK(IWORK),
                LWORK - IWORK + 1, IERR);

            // Copy L to WORK(IU), zeroing out above it

            zlacpy('L', M, M, A, LDA, WORK(IU).asMatrix(LDWRKU), LDWRKU);
            zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero,
                WORK(IU + LDWRKU).asMatrix(LDWRKU), LDWRKU);
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Bidiagonalize L in WORK(IU), copying result to U
            // (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
            // (RWorkspace: need M)

            zgebrd(M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU, S, RWORK(IE),
                WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('L', M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU, U, LDU);

            // Generate right bidiagonalizing vectors in WORK(IU)
            // (CWorkspace: need M*M+3*M, prefer M*M+2*M+(M-1)*NB)
            // (RWorkspace: 0)

            zungbr('P', M, M, M, WORK(IU).asMatrix(LDWRKU), LDWRKU, WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate left bidiagonalizing vectors in U
            // (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB)
            // (RWorkspace: 0)

            zungbr('Q', M, M, M, U, LDU, WORK(ITAUQ), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IRWORK = IE + M;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of L in U and computing right
            // singular vectors of L in WORK(IU)
            // (CWorkspace: need M*M)
            // (RWorkspace: need BDSPAC)

            zbdsqr('U', M, M, M, 0, S, RWORK(IE), WORK(IU).asMatrix(LDWRKU),
                LDWRKU, U, LDU, CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);

            // Multiply right singular vectors of L in WORK(IU) by
            // Q in VT, storing result in A
            // (CWorkspace: need M*M)
            // (RWorkspace: 0)

            zgemm('N', 'N', M, N, M, Complex.one, WORK(IU).asMatrix(LDWRKU),
                LDWRKU, VT, LDVT, Complex.zero, A, LDA);

            // Copy right singular vectors of A from A to VT

            zlacpy('F', M, N, A, LDA, VT, LDVT);
          } else {
            // Insufficient workspace for a fast algorithm

            ITAU = 1;
            IWORK = ITAU + M;

            // Compute A=L*Q, copying result to VT
            // (CWorkspace: need 2*M, prefer M+M*NB)
            // (RWorkspace: 0)

            zgelqf(
                M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR);
            zlacpy('U', M, N, A, LDA, VT, LDVT);

            // Generate Q in VT
            // (CWorkspace: need M+N, prefer M+N*NB)
            // (RWorkspace: 0)

            zunglq(N, N, M, VT, LDVT, WORK(ITAU), WORK(IWORK),
                LWORK - IWORK + 1, IERR);

            // Copy L to U, zeroing out above it

            zlacpy('L', M, M, A, LDA, U, LDU);
            zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero, U(1, 2), LDU);
            IE = 1;
            ITAUQ = ITAU;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Bidiagonalize L in U
            // (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
            // (RWorkspace: need M)

            zgebrd(M, M, U, LDU, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP),
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Multiply right bidiagonalizing vectors in U by Q
            // in VT
            // (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
            // (RWorkspace: 0)

            zunmbr('P', 'L', 'C', M, N, M, U, LDU, WORK(ITAUP), VT, LDVT,
                WORK(IWORK), LWORK - IWORK + 1, IERR);

            // Generate left bidiagonalizing vectors in U
            // (CWorkspace: need 3*M, prefer 2*M+M*NB)
            // (RWorkspace: 0)

            zungbr('Q', M, M, M, U, LDU, WORK(ITAUQ), WORK(IWORK),
                LWORK - IWORK + 1, IERR);
            IRWORK = IE + M;

            // Perform bidiagonal QR iteration, computing left
            // singular vectors of A in U and computing right
            // singular vectors of A in VT
            // (CWorkspace: 0)
            // (RWorkspace: need BDSPAC)

            zbdsqr('U', M, N, M, 0, S, RWORK(IE), VT, LDVT, U, LDU,
                CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
          }
        }
      }
    } else {
      // N < MNTHR

      // Path 10t(N greater than M, but not much larger)
      // Reduce to bidiagonal form without LQ decomposition

      IE = 1;
      ITAUQ = 1;
      ITAUP = ITAUQ + M;
      IWORK = ITAUP + M;

      // Bidiagonalize A
      // (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB)
      // (RWorkspace: M)

      zgebrd(M, N, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP), WORK(IWORK),
          LWORK - IWORK + 1, IERR);
      if (WNTUAS) {
        // If left singular vectors desired in U, copy result to U
        // and generate left bidiagonalizing vectors in U
        // (CWorkspace: need 3*M-1, prefer 2*M+(M-1)*NB)
        // (RWorkspace: 0)

        zlacpy('L', M, M, A, LDA, U, LDU);
        zungbr('Q', M, M, N, U, LDU, WORK(ITAUQ), WORK(IWORK),
            LWORK - IWORK + 1, IERR);
      }
      if (WNTVAS) {
        // If right singular vectors desired in VT, copy result to
        // VT and generate right bidiagonalizing vectors in VT
        // (CWorkspace: need 2*M+NRVT, prefer 2*M+NRVT*NB)
        // (RWorkspace: 0)

        zlacpy('U', M, N, A, LDA, VT, LDVT);
        if (WNTVA) NRVT = N;
        if (WNTVS) NRVT = M;
        zungbr('P', NRVT, N, M, VT, LDVT, WORK(ITAUP), WORK(IWORK),
            LWORK - IWORK + 1, IERR);
      }
      if (WNTUO) {
        // If left singular vectors desired in A, generate left
        // bidiagonalizing vectors in A
        // (CWorkspace: need 3*M-1, prefer 2*M+(M-1)*NB)
        // (RWorkspace: 0)

        zungbr('Q', M, M, N, A, LDA, WORK(ITAUQ), WORK(IWORK),
            LWORK - IWORK + 1, IERR);
      }
      if (WNTVO) {
        // If right singular vectors desired in A, generate right
        // bidiagonalizing vectors in A
        // (CWorkspace: need 3*M, prefer 2*M+M*NB)
        // (RWorkspace: 0)

        zungbr('P', M, N, M, A, LDA, WORK(ITAUP), WORK(IWORK),
            LWORK - IWORK + 1, IERR);
      }
      IRWORK = IE + M;
      if (WNTUAS || WNTUO) NRU = M;
      if (WNTUN) NRU = 0;
      if (WNTVAS || WNTVO) NCVT = N;
      if (WNTVN) NCVT = 0;
      if (!WNTUO && !WNTVO) {
        // Perform bidiagonal QR iteration, if desired, computing
        // left singular vectors in U and computing right singular
        // vectors in VT
        // (CWorkspace: 0)
        // (RWorkspace: need BDSPAC)

        zbdsqr('L', M, NCVT, NRU, 0, S, RWORK(IE), VT, LDVT, U, LDU,
            CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
      } else if (!WNTUO && WNTVO) {
        // Perform bidiagonal QR iteration, if desired, computing
        // left singular vectors in U and computing right singular
        // vectors in A
        // (CWorkspace: 0)
        // (RWorkspace: need BDSPAC)

        zbdsqr('L', M, NCVT, NRU, 0, S, RWORK(IE), A, LDA, U, LDU,
            CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
      } else {
        // Perform bidiagonal QR iteration, if desired, computing
        // left singular vectors in A and computing right singular
        // vectors in VT
        // (CWorkspace: 0)
        // (RWorkspace: need BDSPAC)

        zbdsqr('L', M, NCVT, NRU, 0, S, RWORK(IE), VT, LDVT, A, LDA,
            CDUM.asMatrix(1), 1, RWORK(IRWORK), INFO);
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

  WORK[1] = MAXWRK.toComplex();
}
