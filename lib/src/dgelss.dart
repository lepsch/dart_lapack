import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dbdsqr.dart';
import 'package:lapack/src/dgebrd.dart';
import 'package:lapack/src/dgelqf.dart';
import 'package:lapack/src/dgeqrf.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dorgbr.dart';
import 'package:lapack/src/dormbr.dart';
import 'package:lapack/src/dormlq.dart';
import 'package:lapack/src/dormqr.dart';
import 'package:lapack/src/drscl.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgelss(
  final int M,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Array<double> S_,
  final Box<double> RCOND,
  final Box<int> RANK,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final S = S_.having();
  final WORK = WORK_.having();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  bool LQUERY;
  int BDSPAC,
      BL,
      CHUNK,
      I,
      IASCL,
      IBSCL,
      IE,
      IL,
      ITAU,
      ITAUP,
      ITAUQ,
      IWORK,
      LDWORK,
      MAXMN,
      MAXWRK = 0,
      MINMN,
      MINWRK,
      MM,
      MNTHR = 0;
  int LWORK_DGEQRF,
      LWORK_DORMQR,
      LWORK_DGEBRD,
      LWORK_DORMBR,
      LWORK_DORGBR,
      LWORK_DORMLQ,
      LWORK_DGELQF;
  double ANRM, BIGNUM, BNRM, EPS, SFMIN, SMLNUM, THR;
  final DUM = Array<double>(1);

  // Test the input arguments

  INFO.value = 0;
  MINMN = min(M, N);
  MAXMN = max(M, N);
  LQUERY = (LWORK == -1);
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

  // Compute workspace
  //  (Note: Comments in the code beginning "Workspace:" describe the
  //   minimal amount of workspace needed at that point in the code,
  //   as well as the preferred amount for good performance.
  //   NB refers to the optimal block size for the immediately
  //   following subroutine, as returned by ILAENV.)

  if (INFO.value == 0) {
    MINWRK = 1;
    MAXWRK = 1;
    if (MINMN > 0) {
      MM = M;
      MNTHR = ilaenv(6, 'DGELSS', ' ', M, N, NRHS, -1);
      if (M >= N && M >= MNTHR) {
        // Path 1a - overdetermined, with many more rows than
        //           columns

        // Compute space needed for DGEQRF
        dgeqrf(M, N, A, LDA, DUM(1), DUM(1), -1, INFO);
        LWORK_DGEQRF = DUM[1].toInt();
        // Compute space needed for DORMQR
        dormqr('L', 'T', M, NRHS, N, A, LDA, DUM(1), B, LDB, DUM(1), -1, INFO);
        LWORK_DORMQR = DUM[1].toInt();
        MM = N;
        MAXWRK = max(MAXWRK, N + LWORK_DGEQRF);
        MAXWRK = max(MAXWRK, N + LWORK_DORMQR);
      }
      if (M >= N) {
        // Path 1 - overdetermined or exactly determined

        // Compute workspace needed for DBDSQR

        BDSPAC = max(1, 5 * N);
        // Compute space needed for DGEBRD
        dgebrd(MM, N, A, LDA, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, INFO);
        LWORK_DGEBRD = DUM[1].toInt();
        // Compute space needed for DORMBR
        dormbr('Q', 'L', 'T', MM, NRHS, N, A, LDA, DUM(1), B, LDB, DUM(1), -1,
            INFO);
        LWORK_DORMBR = DUM[1].toInt();
        // Compute space needed for DORGBR
        dorgbr('P', N, N, N, A, LDA, DUM(1), DUM(1), -1, INFO);
        LWORK_DORGBR = DUM[1].toInt();
        // Compute total workspace needed
        MAXWRK = max(MAXWRK, 3 * N + LWORK_DGEBRD);
        MAXWRK = max(MAXWRK, 3 * N + LWORK_DORMBR);
        MAXWRK = max(MAXWRK, 3 * N + LWORK_DORGBR);
        MAXWRK = max(MAXWRK, BDSPAC);
        MAXWRK = max(MAXWRK, N * NRHS);
        MINWRK = max(3 * N + MM, max(3 * N + NRHS, BDSPAC));
        MAXWRK = max(MINWRK, MAXWRK);
      }
      if (N > M) {
        // Compute workspace needed for DBDSQR

        BDSPAC = max(1, 5 * M);
        MINWRK = max(3 * M + NRHS, max(3 * M + N, BDSPAC));
        if (N >= MNTHR) {
          // Path 2a - underdetermined, with many more columns
          // than rows

          // Compute space needed for DGELQF
          dgelqf(M, N, A, LDA, DUM(1), DUM(1), -1, INFO);
          LWORK_DGELQF = DUM[1].toInt();
          // Compute space needed for DGEBRD
          dgebrd(M, M, A, LDA, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, INFO);
          LWORK_DGEBRD = DUM[1].toInt();
          // Compute space needed for DORMBR
          dormbr('Q', 'L', 'T', M, NRHS, N, A, LDA, DUM(1), B, LDB, DUM(1), -1,
              INFO);
          LWORK_DORMBR = DUM[1].toInt();
          // Compute space needed for DORGBR
          dorgbr('P', M, M, M, A, LDA, DUM(1), DUM(1), -1, INFO);
          LWORK_DORGBR = DUM[1].toInt();
          // Compute space needed for DORMLQ
          dormlq(
              'L', 'T', N, NRHS, M, A, LDA, DUM(1), B, LDB, DUM(1), -1, INFO);
          LWORK_DORMLQ = DUM[1].toInt();
          // Compute total workspace needed
          MAXWRK = M + LWORK_DGELQF;
          MAXWRK = max(MAXWRK, M * M + 4 * M + LWORK_DGEBRD);
          MAXWRK = max(MAXWRK, M * M + 4 * M + LWORK_DORMBR);
          MAXWRK = max(MAXWRK, M * M + 4 * M + LWORK_DORGBR);
          MAXWRK = max(MAXWRK, M * M + M + BDSPAC);
          if (NRHS > 1) {
            MAXWRK = max(MAXWRK, M * M + M + M * NRHS);
          } else {
            MAXWRK = max(MAXWRK, M * M + 2 * M);
          }
          MAXWRK = max(MAXWRK, M + LWORK_DORMLQ);
        } else {
          // Path 2 - underdetermined

          // Compute space needed for DGEBRD
          dgebrd(M, N, A, LDA, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, INFO);
          LWORK_DGEBRD = DUM[1].toInt();
          // Compute space needed for DORMBR
          dormbr('Q', 'L', 'T', M, NRHS, M, A, LDA, DUM(1), B, LDB, DUM(1), -1,
              INFO);
          LWORK_DORMBR = DUM[1].toInt();
          // Compute space needed for DORGBR
          dorgbr('P', M, N, M, A, LDA, DUM(1), DUM(1), -1, INFO);
          LWORK_DORGBR = DUM[1].toInt();
          MAXWRK = 3 * M + LWORK_DGEBRD;
          MAXWRK = max(MAXWRK, 3 * M + LWORK_DORMBR);
          MAXWRK = max(MAXWRK, 3 * M + LWORK_DORGBR);
          MAXWRK = max(MAXWRK, BDSPAC);
          MAXWRK = max(MAXWRK, N * NRHS);
        }
      }
      MAXWRK = max(MINWRK, MAXWRK);
    }
    WORK[1] = MAXWRK.toDouble();

    if (LWORK < MINWRK && !LQUERY) INFO.value = -12;
  }

  if (INFO.value != 0) {
    xerbla('DGELSS', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) {
    RANK.value = 0;
    return;
  }

  // Get machine parameters

  EPS = dlamch('P');
  SFMIN = dlamch('S');
  SMLNUM = SFMIN / EPS;
  BIGNUM = ONE / SMLNUM;

  // Scale A if max element outside range [SMLNUM,BIGNUM]

  ANRM = dlange('M', M, N, A, LDA, WORK);
  IASCL = 0;
  if (ANRM > ZERO && ANRM < SMLNUM) {
    // Scale matrix norm up to SMLNUM

    dlascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO);
    IASCL = 1;
  } else if (ANRM > BIGNUM) {
    // Scale matrix norm down to BIGNUM

    dlascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO);
    IASCL = 2;
  } else if (ANRM == ZERO) {
    // Matrix all zero. Return zero solution.

    dlaset('F', max(M, N), NRHS, ZERO, ZERO, B, LDB);
    dlaset('F', MINMN, 1, ZERO, ZERO, S.asMatrix(MINMN), MINMN);
    RANK.value = 0;
    WORK[1] = MAXWRK.toDouble();
    return;
  }

  // Scale B if max element outside range [SMLNUM,BIGNUM]

  BNRM = dlange('M', M, NRHS, B, LDB, WORK);
  IBSCL = 0;
  if (BNRM > ZERO && BNRM < SMLNUM) {
    // Scale matrix norm up to SMLNUM

    dlascl('G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO);
    IBSCL = 1;
  } else if (BNRM > BIGNUM) {
    // Scale matrix norm down to BIGNUM

    dlascl('G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO);
    IBSCL = 2;
  }

  // Overdetermined case

  if (M >= N) {
    // Path 1 - overdetermined or exactly determined

    MM = M;
    if (M >= MNTHR) {
      // Path 1a - overdetermined, with many more rows than columns

      MM = N;
      ITAU = 1;
      IWORK = ITAU + N;

      // Compute A=Q*R
      // (Workspace: need 2*N, prefer N+N*NB)

      dgeqrf(M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, INFO);

      // Multiply B by transpose(Q)
      // (Workspace: need N+NRHS, prefer N+NRHS*NB)

      dormqr('L', 'T', M, NRHS, N, A, LDA, WORK(ITAU), B, LDB, WORK(IWORK),
          LWORK - IWORK + 1, INFO);

      // Zero out below R

      if (N > 1) dlaset('L', N - 1, N - 1, ZERO, ZERO, A(2, 1), LDA);
    }

    IE = 1;
    ITAUQ = IE + N;
    ITAUP = ITAUQ + N;
    IWORK = ITAUP + N;

    // Bidiagonalize R in A
    // (Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB)

    dgebrd(MM, N, A, LDA, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP), WORK(IWORK),
        LWORK - IWORK + 1, INFO);

    // Multiply B by transpose of left bidiagonalizing vectors of R
    // (Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB)

    dormbr('Q', 'L', 'T', MM, NRHS, N, A, LDA, WORK(ITAUQ), B, LDB, WORK(IWORK),
        LWORK - IWORK + 1, INFO);

    // Generate right bidiagonalizing vectors of R in A
    // (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)

    dorgbr('P', N, N, N, A, LDA, WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1,
        INFO);
    IWORK = IE + N;

    // Perform bidiagonal QR iteration
    //   multiply B by transpose of left singular vectors
    //   compute right singular vectors in A
    // (Workspace: need BDSPAC)

    dbdsqr('U', N, N, 0, NRHS, S, WORK(IE), A, LDA, DUM.asMatrix(1), 1, B, LDB,
        WORK(IWORK), INFO);
    if (INFO.value != 0) {
      WORK[1] = MAXWRK.toDouble();
      return;
    }
    // Multiply B by reciprocals of singular values

    THR = max(RCOND.value * S[1], SFMIN);
    if (RCOND.value < ZERO) THR = max(EPS * S[1], SFMIN);
    RANK.value = 0;
    for (I = 1; I <= N; I++) {
      if (S[I] > THR) {
        drscl(NRHS, S[I], B(I, 1).asArray(), LDB);
        RANK.value = RANK.value + 1;
      } else {
        dlaset('F', 1, NRHS, ZERO, ZERO, B(I, 1), LDB);
      }
    }

    // Multiply B by right singular vectors
    // (Workspace: need N, prefer N*NRHS)

    if (LWORK >= LDB * NRHS && NRHS > 1) {
      dgemm('T', 'N', N, NRHS, N, ONE, A, LDA, B, LDB, ZERO, WORK.asMatrix(LDB),
          LDB);
      dlacpy('G', N, NRHS, WORK.asMatrix(LDB), LDB, B, LDB);
    } else if (NRHS > 1) {
      CHUNK = LWORK ~/ N;
      for (I = 1; CHUNK < 0 ? I >= NRHS : I <= NRHS; I += CHUNK) {
        BL = min(NRHS - I + 1, CHUNK);
        dgemm('T', 'N', N, BL, N, ONE, A, LDA, B(1, I), LDB, ZERO,
            WORK.asMatrix(N), N);
        dlacpy('G', N, BL, WORK.asMatrix(N), N, B(1, I), LDB);
      }
    } else if (NRHS == 1) {
      dgemv('T', N, N, ONE, A, LDA, B.asArray(), 1, ZERO, WORK, 1);
      dcopy(N, WORK, 1, B.asArray(), 1);
    }
  } else if (N >= MNTHR &&
      LWORK >= 4 * M + M * M + max(max(M, 2 * M - 4), max(NRHS, N - 3 * M))) {
    // Path 2a - underdetermined, with many more columns than rows
    // and sufficient workspace for an efficient algorithm

    LDWORK = M;
    if (LWORK >=
        max(4 * M + M * LDA + max(max(M, 2 * M - 4), max(NRHS, N - 3 * M)),
            M * LDA + M + M * NRHS)) LDWORK = LDA;
    ITAU = 1;
    IWORK = M + 1;

    // Compute A=L*Q
    // (Workspace: need 2*M, prefer M+M*NB)

    dgelqf(M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, INFO);
    IL = IWORK;

    // Copy L to WORK(IL), zeroing out above it

    dlacpy('L', M, M, A, LDA, WORK(IL).asMatrix(LDWORK), LDWORK);
    dlaset('U', M - 1, M - 1, ZERO, ZERO, WORK(IL + LDWORK).asMatrix(LDWORK),
        LDWORK);
    IE = IL + LDWORK * M;
    ITAUQ = IE + M;
    ITAUP = ITAUQ + M;
    IWORK = ITAUP + M;

    // Bidiagonalize L in WORK(IL)
    // (Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB)

    dgebrd(M, M, WORK(IL).asMatrix(LDWORK), LDWORK, S, WORK(IE), WORK(ITAUQ),
        WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, INFO);

    // Multiply B by transpose of left bidiagonalizing vectors of L
    // (Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)

    dormbr('Q', 'L', 'T', M, NRHS, M, WORK(IL).asMatrix(LDWORK), LDWORK,
        WORK(ITAUQ), B, LDB, WORK(IWORK), LWORK - IWORK + 1, INFO);

    // Generate right bidiagonalizing vectors of R in WORK(IL)
    // (Workspace: need M*M+5*M-1, prefer M*M+4*M+(M-1)*NB)

    dorgbr('P', M, M, M, WORK(IL).asMatrix(LDWORK), LDWORK, WORK(ITAUP),
        WORK(IWORK), LWORK - IWORK + 1, INFO);
    IWORK = IE + M;

    // Perform bidiagonal QR iteration,
    //    computing right singular vectors of L in WORK(IL) and
    //    multiplying B by transpose of left singular vectors
    // (Workspace: need M*M+M+BDSPAC)

    dbdsqr('U', M, M, 0, NRHS, S, WORK(IE), WORK(IL).asMatrix(LDWORK), LDWORK,
        A, LDA, B, LDB, WORK(IWORK), INFO);
    if (INFO.value != 0) {
      WORK[1] = MAXWRK.toDouble();
      return;
    }

    // Multiply B by reciprocals of singular values

    THR = max(RCOND.value * S[1], SFMIN);
    if (RCOND.value < ZERO) THR = max(EPS * S[1], SFMIN);
    RANK.value = 0;
    for (I = 1; I <= M; I++) {
      if (S[I] > THR) {
        drscl(NRHS, S[I], B(I, 1).asArray(), LDB);
        RANK.value = RANK.value + 1;
      } else {
        dlaset('F', 1, NRHS, ZERO, ZERO, B(I, 1), LDB);
      }
    }
    IWORK = IE;

    // Multiply B by right singular vectors of L in WORK(IL)
    // (Workspace: need M*M+2*M, prefer M*M+M+M*NRHS)

    if (LWORK >= LDB * NRHS + IWORK - 1 && NRHS > 1) {
      dgemm('T', 'N', M, NRHS, M, ONE, WORK(IL).asMatrix(LDWORK), LDWORK, B,
          LDB, ZERO, WORK(IWORK).asMatrix(LDB), LDB);
      dlacpy('G', M, NRHS, WORK(IWORK).asMatrix(LDB), LDB, B, LDB);
    } else if (NRHS > 1) {
      CHUNK = (LWORK - IWORK + 1) ~/ M;
      for (I = 1; CHUNK < 0 ? I >= NRHS : I <= NRHS; I += CHUNK) {
        BL = min(NRHS - I + 1, CHUNK);
        dgemm('T', 'N', M, BL, M, ONE, WORK(IL).asMatrix(LDWORK), LDWORK,
            B(1, I), LDB, ZERO, WORK(IWORK).asMatrix(M), M);
        dlacpy('G', M, BL, WORK(IWORK).asMatrix(M), M, B(1, I), LDB);
      }
    } else if (NRHS == 1) {
      dgemv('T', M, M, ONE, WORK(IL).asMatrix(LDWORK), LDWORK,
          B(1, 1).asArray(), 1, ZERO, WORK(IWORK), 1);
      dcopy(M, WORK(IWORK), 1, B(1, 1).asArray(), 1);
    }

    // Zero out below first M rows of B

    dlaset('F', N - M, NRHS, ZERO, ZERO, B(M + 1, 1), LDB);
    IWORK = ITAU + M;

    // Multiply transpose(Q) by B
    // (Workspace: need M+NRHS, prefer M+NRHS*NB)

    dormlq('L', 'T', N, NRHS, M, A, LDA, WORK(ITAU), B, LDB, WORK(IWORK),
        LWORK - IWORK + 1, INFO);
  } else {
    // Path 2 - remaining underdetermined cases

    IE = 1;
    ITAUQ = IE + M;
    ITAUP = ITAUQ + M;
    IWORK = ITAUP + M;

    // Bidiagonalize A
    // (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)

    dgebrd(M, N, A, LDA, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP), WORK(IWORK),
        LWORK - IWORK + 1, INFO);

    // Multiply B by transpose of left bidiagonalizing vectors
    // (Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB)

    dormbr('Q', 'L', 'T', M, NRHS, N, A, LDA, WORK(ITAUQ), B, LDB, WORK(IWORK),
        LWORK - IWORK + 1, INFO);

    // Generate right bidiagonalizing vectors in A
    // (Workspace: need 4*M, prefer 3*M+M*NB)

    dorgbr('P', M, N, M, A, LDA, WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1,
        INFO);
    IWORK = IE + M;

    // Perform bidiagonal QR iteration,
    //    computing right singular vectors of A in A and
    //    multiplying B by transpose of left singular vectors
    // (Workspace: need BDSPAC)

    dbdsqr('L', M, N, 0, NRHS, S, WORK(IE), A, LDA, DUM.asMatrix(1), 1, B, LDB,
        WORK(IWORK), INFO);
    if (INFO.value != 0) {
      WORK[1] = MAXWRK.toDouble();
      return;
    }

    // Multiply B by reciprocals of singular values

    THR = max(RCOND.value * S[1], SFMIN);
    if (RCOND.value < ZERO) THR = max(EPS * S[1], SFMIN);
    RANK.value = 0;
    for (I = 1; I <= M; I++) {
      if (S[I] > THR) {
        drscl(NRHS, S[I], B(I, 1).asArray(), LDB);
        RANK.value = RANK.value + 1;
      } else {
        dlaset('F', 1, NRHS, ZERO, ZERO, B(I, 1), LDB);
      }
    }

    // Multiply B by right singular vectors of A
    // (Workspace: need N, prefer N*NRHS)

    if (LWORK >= LDB * NRHS && NRHS > 1) {
      dgemm('T', 'N', N, NRHS, M, ONE, A, LDA, B, LDB, ZERO, WORK.asMatrix(LDB),
          LDB);
      dlacpy('F', N, NRHS, WORK.asMatrix(LDB), LDB, B, LDB);
    } else if (NRHS > 1) {
      CHUNK = LWORK ~/ N;
      for (I = 1; CHUNK < 0 ? I >= NRHS : I <= NRHS; I += CHUNK) {
        BL = min(NRHS - I + 1, CHUNK);
        dgemm('T', 'N', N, BL, M, ONE, A, LDA, B(1, I), LDB, ZERO,
            WORK.asMatrix(N), N);
        dlacpy('F', N, BL, WORK.asMatrix(N), N, B(1, I), LDB);
      }
    } else if (NRHS == 1) {
      dgemv('T', M, N, ONE, A, LDA, B.asArray(), 1, ZERO, WORK, 1);
      dcopy(N, WORK, 1, B.asArray(), 1);
    }
  }

  // Undo scaling

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
}
