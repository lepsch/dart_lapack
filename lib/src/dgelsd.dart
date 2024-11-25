// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/lapack.dart';

void dgelsd(
  final int M,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Array<double> S_,
  final double RCOND,
  final Box<int> RANK,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final S = S_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;

  // Test the input arguments.
  INFO.value = 0;
  final MINMN = max(1, min(M, N));
  final MAXMN = max(M, N);
  final MNTHR = ilaenv(6, 'DGELSD', ' ', M, N, NRHS, -1);
  final LQUERY = LWORK == -1;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, M)) {
    INFO.value = -5;
  } else if (LDB < max(1, MAXMN)) {
    INFO.value = -7;
  }

  if (INFO.value != 0) {
    xerbla('DGELSD', -INFO.value);
    return;
  }

  final SMLSIZ = ilaenv(9, 'DGELSD', ' ', 0, 0, 0, 0);

  // Compute workspace.
  // (Note: Comments in the code beginning "Workspace:" describe the
  // minimal amount of workspace needed at that point in the code,
  // as well as the preferred amount for good performance.
  // NB refers to the optimal block size for the immediately
  // following subroutine, as returned by ILAENV.)

  final NLVL = max(log(MINMN / (SMLSIZ + 1)) ~/ log(TWO) + 1, 0);
  final LIWORK = 3 * MINMN * NLVL + 11 * MINMN;
  int MAXWRK = 1;
  final int MM;
  if (M >= N && M >= MNTHR) {
    // Path 1a - overdetermined, with many more rows than columns.
    MM = N;
    MAXWRK = max(MAXWRK, N + N * ilaenv(1, 'DGEQRF', ' ', M, N, -1, -1));
    MAXWRK = max(MAXWRK, N + NRHS * ilaenv(1, 'DORMQR', 'LT', M, NRHS, N, -1));
  } else {
    MM = M;
  }
  int MINWRK = 1;
  final int WLALSD;
  if (M >= N) {
    // Path 1 - overdetermined or exactly determined.
    MAXWRK =
        max(MAXWRK, 3 * N + (MM + N) * ilaenv(1, 'DGEBRD', ' ', MM, N, -1, -1));
    MAXWRK =
        max(MAXWRK, 3 * N + NRHS * ilaenv(1, 'DORMBR', 'QLT', MM, NRHS, N, -1));
    MAXWRK = max(
        MAXWRK, 3 * N + (N - 1) * ilaenv(1, 'DORMBR', 'PLN', N, NRHS, N, -1));
    WLALSD = 9 * N +
        2 * N * SMLSIZ +
        8 * N * NLVL +
        N * NRHS +
        pow((SMLSIZ + 1), 2).toInt();
    MAXWRK = max(MAXWRK, 3 * N + WLALSD);
    MINWRK = max(3 * N + MM, max(3 * N + NRHS, 3 * N + WLALSD));
  } else {
    WLALSD = 9 * M +
        2 * M * SMLSIZ +
        8 * M * NLVL +
        M * NRHS +
        pow((SMLSIZ + 1), 2).toInt();
    if (N >= MNTHR) {
      // Path 2a - underdetermined, with many more columns than rows.
      MAXWRK = M + M * ilaenv(1, 'DGELQF', ' ', M, N, -1, -1);
      MAXWRK = max(MAXWRK,
          M * M + 4 * M + 2 * M * ilaenv(1, 'DGEBRD', ' ', M, M, -1, -1));
      MAXWRK = max(MAXWRK,
          M * M + 4 * M + NRHS * ilaenv(1, 'DORMBR', 'QLT', M, NRHS, M, -1));
      MAXWRK = max(MAXWRK,
          M * M + 4 * M + (M - 1) * ilaenv(1, 'DORMBR', 'PLN', M, NRHS, M, -1));
      if (NRHS > 1) {
        MAXWRK = max(MAXWRK, M * M + M + M * NRHS);
      } else {
        MAXWRK = max(MAXWRK, M * M + 2 * M);
      }
      MAXWRK =
          max(MAXWRK, M + NRHS * ilaenv(1, 'DORMLQ', 'LT', N, NRHS, M, -1));
      MAXWRK = max(MAXWRK, M * M + 4 * M + WLALSD);
      // XXX: Ensure the Path 2a case below is triggered.  The workspace
      // calculation should use queries for all routines eventually.
      MAXWRK = max(
          MAXWRK, 4 * M + M * M + max(max(M, 2 * M - 4), max(NRHS, N - 3 * M)));
    } else {
      // Path 2 - remaining underdetermined cases.
      MAXWRK = 3 * M + (N + M) * ilaenv(1, 'DGEBRD', ' ', M, N, -1, -1);
      MAXWRK = max(
          MAXWRK, 3 * M + NRHS * ilaenv(1, 'DORMBR', 'QLT', M, NRHS, N, -1));
      MAXWRK =
          max(MAXWRK, 3 * M + M * ilaenv(1, 'DORMBR', 'PLN', N, NRHS, M, -1));
      MAXWRK = max(MAXWRK, 3 * M + WLALSD);
    }
    MINWRK = max(3 * M + NRHS, max(3 * M + M, 3 * M + WLALSD));
  }
  MINWRK = min(MINWRK, MAXWRK);
  WORK[1] = MAXWRK.toDouble();
  IWORK[1] = LIWORK;

  if (LWORK < MINWRK && !LQUERY) {
    INFO.value = -12;
    xerbla('DGELSD', -INFO.value);
    return;
  } else if (LQUERY) {
    WORK[1] = MAXWRK.toDouble();
    IWORK[1] = LIWORK;
    return;
  }

  // Quick return if possible.
  if (M == 0 || N == 0) {
    RANK.value = 0;
    return;
  }

  // Get machine parameters.
  final EPS = dlamch('P');
  final SFMIN = dlamch('S');
  final SMLNUM = SFMIN / EPS;
  final BIGNUM = ONE / SMLNUM;

  // Scale A if max entry outside range [SMLNUM,BIGNUM].
  final ANRM = dlange('M', M, N, A, LDA, WORK);
  final int IASCL;
  if (ANRM > ZERO && ANRM < SMLNUM) {
    // Scale matrix norm up to SMLNUM.

    dlascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO);
    IASCL = 1;
  } else if (ANRM > BIGNUM) {
    // Scale matrix norm down to BIGNUM.

    dlascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO);
    IASCL = 2;
  } else if (ANRM == ZERO) {
    // Matrix all zero. Return zero solution.

    dlaset('F', max(M, N), NRHS, ZERO, ZERO, B, LDB);
    dlaset('F', MINMN, 1, ZERO, ZERO, S.asMatrix(1), 1);
    RANK.value = 0;

    WORK[1] = MAXWRK.toDouble();
    IWORK[1] = LIWORK;
    return;
  } else {
    IASCL = 0;
  }

  // Scale B if max entry outside range [SMLNUM,BIGNUM].
  final BNRM = dlange('M', M, NRHS, B, LDB, WORK);
  final int IBSCL;
  if (BNRM > ZERO && BNRM < SMLNUM) {
    // Scale matrix norm up to SMLNUM.
    dlascl('G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO);
    IBSCL = 1;
  } else if (BNRM > BIGNUM) {
    // Scale matrix norm down to BIGNUM.
    dlascl('G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO);
    IBSCL = 2;
  } else {
    IBSCL = 0;
  }

  // If M < N make sure certain entries of B are zero.
  if (M < N) dlaset('F', N - M, NRHS, ZERO, ZERO, B(M + 1, 1), LDB);

  // Overdetermined case.
  if (M >= N) {
    // Path 1 - overdetermined or exactly determined.
    final int MM;
    if (M >= MNTHR) {
      // Path 1a - overdetermined, with many more rows than columns.
      MM = N;
      final ITAU = 1;
      final NWORK = ITAU + N;

      // Compute A=Q*R.
      // (Workspace: need 2*N, prefer N+N*NB)
      dgeqrf(M, N, A, LDA, WORK(ITAU), WORK(NWORK), LWORK - NWORK + 1, INFO);

      // Multiply B by transpose(Q).
      // (Workspace: need N+NRHS, prefer N+NRHS*NB)
      dormqr('L', 'T', M, NRHS, N, A, LDA, WORK(ITAU), B, LDB, WORK(NWORK),
          LWORK - NWORK + 1, INFO);

      // Zero out below R.
      if (N > 1) {
        dlaset('L', N - 1, N - 1, ZERO, ZERO, A(2, 1), LDA);
      }
    } else {
      MM = M;
    }

    final IE = 1;
    final ITAUQ = IE + N;
    final ITAUP = ITAUQ + N;
    final NWORK = ITAUP + N;

    // Bidiagonalize R in A.
    // (Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB)
    dgebrd(MM, N, A, LDA, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP), WORK(NWORK),
        LWORK - NWORK + 1, INFO);

    // Multiply B by transpose of left bidiagonalizing vectors of R.
    // (Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB)
    dormbr('Q', 'L', 'T', MM, NRHS, N, A, LDA, WORK(ITAUQ), B, LDB, WORK(NWORK),
        LWORK - NWORK + 1, INFO);

    // Solve the bidiagonal least squares problem.
    dlalsd('U', SMLSIZ, N, NRHS, S, WORK(IE), B, LDB, RCOND, RANK, WORK(NWORK),
        IWORK, INFO);
    if (INFO.value != 0) {
      WORK[1] = MAXWRK.toDouble();
      IWORK[1] = LIWORK;
      return;
    }

    // Multiply B by right bidiagonalizing vectors of R.
    dormbr('P', 'L', 'N', N, NRHS, N, A, LDA, WORK(ITAUP), B, LDB, WORK(NWORK),
        LWORK - NWORK + 1, INFO);
  } else if (N >= MNTHR &&
      LWORK >=
          4 * M +
              M * M +
              max(max(M, 2 * M - 4), max(max(NRHS, N - 3 * M), WLALSD))) {
    // Path 2a - underdetermined, with many more columns than rows
    // and sufficient workspace for an efficient algorithm.
    final LDWORK = (LWORK >=
            max(4 * M + M * LDA + max(max(M, 2 * M - 4), max(NRHS, N - 3 * M)),
                max(M * LDA + M + M * NRHS, 4 * M + M * LDA + WLALSD)))
        ? LDA
        : M;
    final ITAU = 1;
    var NWORK = M + 1;

    // Compute A=L*Q.
    // (Workspace: need 2*M, prefer M+M*NB)
    dgelqf(M, N, A, LDA, WORK(ITAU), WORK(NWORK), LWORK - NWORK + 1, INFO);
    final IL = NWORK;

    // Copy L to WORK(IL), zeroing out above its diagonal.
    dlacpy('L', M, M, A, LDA, WORK(IL).asMatrix(LDWORK), LDWORK);
    dlaset('U', M - 1, M - 1, ZERO, ZERO, WORK(IL + LDWORK).asMatrix(LDWORK),
        LDWORK);
    final IE = IL + LDWORK * M;
    final ITAUQ = IE + M;
    final ITAUP = ITAUQ + M;
    NWORK = ITAUP + M;

    // Bidiagonalize L in WORK(IL).
    // (Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB)
    dgebrd(M, M, WORK(IL).asMatrix(LDWORK), LDWORK, S, WORK(IE), WORK(ITAUQ),
        WORK(ITAUP), WORK(NWORK), LWORK - NWORK + 1, INFO);

    // Multiply B by transpose of left bidiagonalizing vectors of L.
    // (Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)
    dormbr('Q', 'L', 'T', M, NRHS, M, WORK(IL).asMatrix(LDWORK), LDWORK,
        WORK(ITAUQ), B, LDB, WORK(NWORK), LWORK - NWORK + 1, INFO);

    // Solve the bidiagonal least squares problem.
    dlalsd('U', SMLSIZ, M, NRHS, S, WORK(IE), B, LDB, RCOND, RANK, WORK(NWORK),
        IWORK, INFO);
    if (INFO.value != 0) {
      WORK[1] = MAXWRK.toDouble();
      IWORK[1] = LIWORK;
      return;
    }

    // Multiply B by right bidiagonalizing vectors of L.
    dormbr('P', 'L', 'N', M, NRHS, M, WORK(IL).asMatrix(LDWORK), LDWORK,
        WORK(ITAUP), B, LDB, WORK(NWORK), LWORK - NWORK + 1, INFO);

    // Zero out below first M rows of B.
    dlaset('F', N - M, NRHS, ZERO, ZERO, B(M + 1, 1), LDB);
    NWORK = ITAU + M;

    // Multiply transpose(Q) by B.
    // (Workspace: need M+NRHS, prefer M+NRHS*NB)
    dormlq('L', 'T', N, NRHS, M, A, LDA, WORK(ITAU), B, LDB, WORK(NWORK),
        LWORK - NWORK + 1, INFO);
  } else {
    // Path 2 - remaining underdetermined cases.
    final IE = 1;
    final ITAUQ = IE + M;
    final ITAUP = ITAUQ + M;
    final NWORK = ITAUP + M;

    // Bidiagonalize A.
    // (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)
    dgebrd(M, N, A, LDA, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP), WORK(NWORK),
        LWORK - NWORK + 1, INFO);

    // Multiply B by transpose of left bidiagonalizing vectors.
    // (Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB)
    dormbr('Q', 'L', 'T', M, NRHS, N, A, LDA, WORK(ITAUQ), B, LDB, WORK(NWORK),
        LWORK - NWORK + 1, INFO);

    // Solve the bidiagonal least squares problem.
    dlalsd('L', SMLSIZ, M, NRHS, S, WORK(IE), B, LDB, RCOND, RANK, WORK(NWORK),
        IWORK, INFO);
    if (INFO.value != 0) {
      WORK[1] = MAXWRK.toDouble();
      IWORK[1] = LIWORK;
      return;
    }

    // Multiply B by right bidiagonalizing vectors of A.
    dormbr('P', 'L', 'N', N, NRHS, M, A, LDA, WORK(ITAUP), B, LDB, WORK(NWORK),
        LWORK - NWORK + 1, INFO);
  }

  // Undo scaling.
  if (IASCL == 1) {
    dlascl('G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO);
    dlascl('G', 0, 0, SMLNUM, ANRM, MINMN, 1, S.asMatrix(MINMN), MINMN, INFO);
  } else if (IASCL == 2) {
    dlascl('G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO);
    dlascl('G', 0, 0, BIGNUM, ANRM, MINMN, 1, S.asMatrix(MINMN), MINMN, INFO);
  }
  if (IBSCL == 1) {
    dlascl('G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO);
  } else if (IBSCL == 2) {
    dlascl('G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO);
  }
  WORK[1] = MAXWRK.toDouble();
  IWORK[1] = LIWORK;
}
