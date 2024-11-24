// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zgeqrf.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zbdsqr.dart';
import 'package:lapack/src/zdrscl.dart';
import 'package:lapack/src/zgebrd.dart';
import 'package:lapack/src/zgelqf.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlascl.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zungbr.dart';
import 'package:lapack/src/zunmbr.dart';
import 'package:lapack/src/zunmlq.dart';
import 'package:lapack/src/zunmqr.dart';

void zgelss(
  final int M,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<double> S_,
  final double RCOND,
  final Box<int> RANK,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final S = S_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool LQUERY;
  int BL,
      CHUNK,
      I,
      IASCL,
      IBSCL,
      IE,
      IL,
      IRWORK,
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
  int LWORK_ZGEBRD, LWORK_ZUNMBR, LWORK_ZUNGBR, LWORK_ZUNMLQ, LWORK_ZGELQF;
  double ANRM, BIGNUM, BNRM, EPS, SFMIN, SMLNUM, THR;
  final DUM = Array<Complex>(1);

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
  //   CWorkspace refers to complex workspace, and RWorkspace refers
  //   to real workspace. NB refers to the optimal block size for the
  //   immediately following subroutine, as returned by ILAENV.)

  if (INFO.value == 0) {
    MINWRK = 1;
    MAXWRK = 1;
    if (MINMN > 0) {
      MM = M;
      MNTHR = ilaenv(6, 'ZGELSS', ' ', M, N, NRHS, -1);
      if (M >= N && M >= MNTHR) {
        // Path 1a - overdetermined, with many more rows than
        //           columns

        // Compute space needed for ZGEQRF
        zgeqrf(M, N, A, LDA, DUM(1), DUM(1), -1, INFO);
        // Compute space needed for ZUNMQR
        zunmqr('L', 'C', M, NRHS, N, A, LDA, DUM(1), B, LDB, DUM(1), -1, INFO);
        MM = N;
        MAXWRK = max(MAXWRK, N + N * ilaenv(1, 'ZGEQRF', ' ', M, N, -1, -1));
        MAXWRK =
            max(MAXWRK, N + NRHS * ilaenv(1, 'ZUNMQR', 'LC', M, NRHS, N, -1));
      }
      if (M >= N) {
        // Path 1 - overdetermined or exactly determined

        // Compute space needed for ZGEBRD
        zgebrd(MM, N, A, LDA, S, S, DUM(1), DUM(1), DUM(1), -1, INFO);
        LWORK_ZGEBRD = DUM[1].toInt();
        // Compute space needed for ZUNMBR
        zunmbr('Q', 'L', 'C', MM, NRHS, N, A, LDA, DUM(1), B, LDB, DUM(1), -1,
            INFO);
        LWORK_ZUNMBR = DUM[1].toInt();
        // Compute space needed for ZUNGBR
        zungbr('P', N, N, N, A, LDA, DUM(1), DUM(1), -1, INFO);
        LWORK_ZUNGBR = DUM[1].toInt();
        // Compute total workspace needed
        MAXWRK = max(MAXWRK, 2 * N + LWORK_ZGEBRD);
        MAXWRK = max(MAXWRK, 2 * N + LWORK_ZUNMBR);
        MAXWRK = max(MAXWRK, 2 * N + LWORK_ZUNGBR);
        MAXWRK = max(MAXWRK, N * NRHS);
        MINWRK = 2 * N + max(NRHS, M);
      }
      if (N > M) {
        MINWRK = 2 * M + max(NRHS, N);
        if (N >= MNTHR) {
          // Path 2a - underdetermined, with many more columns
          // than rows

          // Compute space needed for ZGELQF
          zgelqf(M, N, A, LDA, DUM(1), DUM(1), -1, INFO);
          LWORK_ZGELQF = DUM[1].toInt();
          // Compute space needed for ZGEBRD
          zgebrd(M, M, A, LDA, S, S, DUM(1), DUM(1), DUM(1), -1, INFO);
          LWORK_ZGEBRD = DUM[1].toInt();
          // Compute space needed for ZUNMBR
          zunmbr('Q', 'L', 'C', M, NRHS, N, A, LDA, DUM(1), B, LDB, DUM(1), -1,
              INFO);
          LWORK_ZUNMBR = DUM[1].toInt();
          // Compute space needed for ZUNGBR
          zungbr('P', M, M, M, A, LDA, DUM(1), DUM(1), -1, INFO);
          LWORK_ZUNGBR = DUM[1].toInt();
          // Compute space needed for ZUNMLQ
          zunmlq(
              'L', 'C', N, NRHS, M, A, LDA, DUM(1), B, LDB, DUM(1), -1, INFO);
          LWORK_ZUNMLQ = DUM[1].toInt();
          // Compute total workspace needed
          MAXWRK = M + LWORK_ZGELQF;
          MAXWRK = max(MAXWRK, 3 * M + M * M + LWORK_ZGEBRD);
          MAXWRK = max(MAXWRK, 3 * M + M * M + LWORK_ZUNMBR);
          MAXWRK = max(MAXWRK, 3 * M + M * M + LWORK_ZUNGBR);
          if (NRHS > 1) {
            MAXWRK = max(MAXWRK, M * M + M + M * NRHS);
          } else {
            MAXWRK = max(MAXWRK, M * M + 2 * M);
          }
          MAXWRK = max(MAXWRK, M + LWORK_ZUNMLQ);
        } else {
          // Path 2 - underdetermined

          // Compute space needed for ZGEBRD
          zgebrd(M, N, A, LDA, S, S, DUM(1), DUM(1), DUM(1), -1, INFO);
          LWORK_ZGEBRD = DUM[1].toInt();
          // Compute space needed for ZUNMBR
          zunmbr('Q', 'L', 'C', M, NRHS, M, A, LDA, DUM(1), B, LDB, DUM(1), -1,
              INFO);
          LWORK_ZUNMBR = DUM[1].toInt();
          // Compute space needed for ZUNGBR
          zungbr('P', M, N, M, A, LDA, DUM(1), DUM(1), -1, INFO);
          LWORK_ZUNGBR = DUM[1].toInt();
          MAXWRK = 2 * M + LWORK_ZGEBRD;
          MAXWRK = max(MAXWRK, 2 * M + LWORK_ZUNMBR);
          MAXWRK = max(MAXWRK, 2 * M + LWORK_ZUNGBR);
          MAXWRK = max(MAXWRK, N * NRHS);
        }
      }
      MAXWRK = max(MINWRK, MAXWRK);
    }
    WORK[1] = MAXWRK.toComplex();

    if (LWORK < MINWRK && !LQUERY) INFO.value = -12;
  }

  if (INFO.value != 0) {
    xerbla('ZGELSS', -INFO.value);
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

  ANRM = zlange('M', M, N, A, LDA, RWORK);
  IASCL = 0;
  if (ANRM > ZERO && ANRM < SMLNUM) {
    // Scale matrix norm up to SMLNUM

    zlascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO);
    IASCL = 1;
  } else if (ANRM > BIGNUM) {
    // Scale matrix norm down to BIGNUM

    zlascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO);
    IASCL = 2;
  } else if (ANRM == ZERO) {
    // Matrix all zero. Return zero solution.

    zlaset('F', max(M, N), NRHS, Complex.zero, Complex.zero, B, LDB);
    dlaset('F', MINMN, 1, ZERO, ZERO, S.asMatrix(MINMN), MINMN);
    RANK.value = 0;
    WORK[1] = MAXWRK.toComplex();
    return;
  }

  // Scale B if max element outside range [SMLNUM,BIGNUM]

  BNRM = zlange('M', M, NRHS, B, LDB, RWORK);
  IBSCL = 0;
  if (BNRM > ZERO && BNRM < SMLNUM) {
    // Scale matrix norm up to SMLNUM

    zlascl('G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO);
    IBSCL = 1;
  } else if (BNRM > BIGNUM) {
    // Scale matrix norm down to BIGNUM

    zlascl('G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO);
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
      // (CWorkspace: need 2*N, prefer N+N*NB)
      // (RWorkspace: none)

      zgeqrf(M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, INFO);

      // Multiply B by transpose(Q)
      // (CWorkspace: need N+NRHS, prefer N+NRHS*NB)
      // (RWorkspace: none)

      zunmqr('L', 'C', M, NRHS, N, A, LDA, WORK(ITAU), B, LDB, WORK(IWORK),
          LWORK - IWORK + 1, INFO);

      // Zero out below R

      if (N > 1) {
        zlaset('L', N - 1, N - 1, Complex.zero, Complex.zero, A(2, 1), LDA);
      }
    }

    IE = 1;
    ITAUQ = 1;
    ITAUP = ITAUQ + N;
    IWORK = ITAUP + N;

    // Bidiagonalize R in A
    // (CWorkspace: need 2*N+MM, prefer 2*N+(MM+N)*NB)
    // (RWorkspace: need N)

    zgebrd(MM, N, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP), WORK(IWORK),
        LWORK - IWORK + 1, INFO);

    // Multiply B by transpose of left bidiagonalizing vectors of R
    // (CWorkspace: need 2*N+NRHS, prefer 2*N+NRHS*NB)
    // (RWorkspace: none)

    zunmbr('Q', 'L', 'C', MM, NRHS, N, A, LDA, WORK(ITAUQ), B, LDB, WORK(IWORK),
        LWORK - IWORK + 1, INFO);

    // Generate right bidiagonalizing vectors of R in A
    // (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
    // (RWorkspace: none)

    zungbr('P', N, N, N, A, LDA, WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1,
        INFO);
    IRWORK = IE + N;

    // Perform bidiagonal QR iteration
    //   multiply B by transpose of left singular vectors
    //   compute right singular vectors in A
    // (CWorkspace: none)
    // (RWorkspace: need BDSPAC)

    zbdsqr('U', N, N, 0, NRHS, S, RWORK(IE), A, LDA, DUM.asMatrix(1), 1, B, LDB,
        RWORK(IRWORK), INFO);
    if (INFO.value != 0) {
      WORK[1] = MAXWRK.toComplex();
      return;
    }

    // Multiply B by reciprocals of singular values

    THR = max(RCOND * S[1], SFMIN);
    if (RCOND < ZERO) THR = max(EPS * S[1], SFMIN);
    RANK.value = 0;
    for (I = 1; I <= N; I++) {
      if (S[I] > THR) {
        zdrscl(NRHS, S[I], B(I, 1).asArray(), LDB);
        RANK.value++;
      } else {
        zlaset('F', 1, NRHS, Complex.zero, Complex.zero, B(I, 1), LDB);
      }
    }

    // Multiply B by right singular vectors
    // (CWorkspace: need N, prefer N*NRHS)
    // (RWorkspace: none)

    if (LWORK >= LDB * NRHS && NRHS > 1) {
      zgemm('C', 'N', N, NRHS, N, Complex.one, A, LDA, B, LDB, Complex.zero,
          WORK.asMatrix(LDB), LDB);
      zlacpy('G', N, NRHS, WORK.asMatrix(LDB), LDB, B, LDB);
    } else if (NRHS > 1) {
      CHUNK = LWORK ~/ N;
      for (I = 1; I <= NRHS; I += CHUNK) {
        BL = min(NRHS - I + 1, CHUNK);
        zgemm('C', 'N', N, BL, N, Complex.one, A, LDA, B(1, I), LDB,
            Complex.zero, WORK.asMatrix(N), N);
        zlacpy('G', N, BL, WORK.asMatrix(N), N, B(1, I), LDB);
      }
    } else if (NRHS == 1) {
      zgemv('C', N, N, Complex.one, A, LDA, B.asArray(), 1, Complex.zero, WORK,
          1);
      zcopy(N, WORK, 1, B.asArray(), 1);
    }
  } else if (N >= MNTHR &&
      LWORK >= 3 * M + M * M + max(M, max(NRHS, N - 2 * M))) {
    // Underdetermined case, M much less than N

    // Path 2a - underdetermined, with many more columns than rows
    // and sufficient workspace for an efficient algorithm

    LDWORK = M;
    if (LWORK >= 3 * M + M * LDA + max(M, max(NRHS, N - 2 * M))) LDWORK = LDA;
    ITAU = 1;
    IWORK = M + 1;

    // Compute A=L*Q
    // (CWorkspace: need 2*M, prefer M+M*NB)
    // (RWorkspace: none)

    zgelqf(M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, INFO);
    IL = IWORK;

    // Copy L to WORK(IL), zeroing out above it

    zlacpy('L', M, M, A, LDA, WORK(IL).asMatrix(LDWORK), LDWORK);
    zlaset('U', M - 1, M - 1, Complex.zero, Complex.zero,
        WORK(IL + LDWORK).asMatrix(LDWORK), LDWORK);
    IE = 1;
    ITAUQ = IL + LDWORK * M;
    ITAUP = ITAUQ + M;
    IWORK = ITAUP + M;

    // Bidiagonalize L in WORK(IL)
    // (CWorkspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
    // (RWorkspace: need M)

    zgebrd(M, M, WORK(IL).asMatrix(LDWORK), LDWORK, S, RWORK(IE), WORK(ITAUQ),
        WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, INFO);

    // Multiply B by transpose of left bidiagonalizing vectors of L
    // (CWorkspace: need M*M+3*M+NRHS, prefer M*M+3*M+NRHS*NB)
    // (RWorkspace: none)

    zunmbr('Q', 'L', 'C', M, NRHS, M, WORK(IL).asMatrix(LDWORK), LDWORK,
        WORK(ITAUQ), B, LDB, WORK(IWORK), LWORK - IWORK + 1, INFO);

    // Generate right bidiagonalizing vectors of R in WORK(IL)
    // (CWorkspace: need M*M+4*M-1, prefer M*M+3*M+(M-1)*NB)
    // (RWorkspace: none)

    zungbr('P', M, M, M, WORK(IL).asMatrix(LDWORK), LDWORK, WORK(ITAUP),
        WORK(IWORK), LWORK - IWORK + 1, INFO);
    IRWORK = IE + M;

    // Perform bidiagonal QR iteration, computing right singular
    // vectors of L in WORK(IL) and multiplying B by transpose of
    // left singular vectors
    // (CWorkspace: need M*M)
    // (RWorkspace: need BDSPAC)

    zbdsqr('U', M, M, 0, NRHS, S, RWORK(IE), WORK(IL).asMatrix(LDWORK), LDWORK,
        A, LDA, B, LDB, RWORK(IRWORK), INFO);
    if (INFO.value != 0) {
      WORK[1] = MAXWRK.toComplex();
      return;
    }

    // Multiply B by reciprocals of singular values

    THR = max(RCOND * S[1], SFMIN);
    if (RCOND < ZERO) THR = max(EPS * S[1], SFMIN);
    RANK.value = 0;
    for (I = 1; I <= M; I++) {
      if (S[I] > THR) {
        zdrscl(NRHS, S[I], B(I, 1).asArray(), LDB);
        RANK.value++;
      } else {
        zlaset('F', 1, NRHS, Complex.zero, Complex.zero, B(I, 1), LDB);
      }
    }
    IWORK = IL + M * LDWORK;

    // Multiply B by right singular vectors of L in WORK(IL)
    // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NRHS)
    // (RWorkspace: none)

    if (LWORK >= LDB * NRHS + IWORK - 1 && NRHS > 1) {
      zgemm('C', 'N', M, NRHS, M, Complex.one, WORK(IL).asMatrix(LDWORK),
          LDWORK, B, LDB, Complex.zero, WORK(IWORK).asMatrix(LDB), LDB);
      zlacpy('G', M, NRHS, WORK(IWORK).asMatrix(LDB), LDB, B, LDB);
    } else if (NRHS > 1) {
      CHUNK = (LWORK - IWORK + 1) ~/ M;
      for (I = 1; I <= NRHS; I += CHUNK) {
        BL = min(NRHS - I + 1, CHUNK);
        zgemm('C', 'N', M, BL, M, Complex.one, WORK(IL).asMatrix(LDWORK),
            LDWORK, B(1, I), LDB, Complex.zero, WORK(IWORK).asMatrix(M), M);
        zlacpy('G', M, BL, WORK(IWORK).asMatrix(M), M, B(1, I), LDB);
      }
    } else if (NRHS == 1) {
      zgemv('C', M, M, Complex.one, WORK(IL).asMatrix(LDWORK), LDWORK,
          B(1, 1).asArray(), 1, Complex.zero, WORK(IWORK), 1);
      zcopy(M, WORK(IWORK), 1, B(1, 1).asArray(), 1);
    }

    // Zero out below first M rows of B

    zlaset('F', N - M, NRHS, Complex.zero, Complex.zero, B(M + 1, 1), LDB);
    IWORK = ITAU + M;

    // Multiply transpose(Q) by B
    // (CWorkspace: need M+NRHS, prefer M+NHRS*NB)
    // (RWorkspace: none)

    zunmlq('L', 'C', N, NRHS, M, A, LDA, WORK(ITAU), B, LDB, WORK(IWORK),
        LWORK - IWORK + 1, INFO);
  } else {
    // Path 2 - remaining underdetermined cases

    IE = 1;
    ITAUQ = 1;
    ITAUP = ITAUQ + M;
    IWORK = ITAUP + M;

    // Bidiagonalize A
    // (CWorkspace: need 3*M, prefer 2*M+(M+N)*NB)
    // (RWorkspace: need N)

    zgebrd(M, N, A, LDA, S, RWORK(IE), WORK(ITAUQ), WORK(ITAUP), WORK(IWORK),
        LWORK - IWORK + 1, INFO);

    // Multiply B by transpose of left bidiagonalizing vectors
    // (CWorkspace: need 2*M+NRHS, prefer 2*M+NRHS*NB)
    // (RWorkspace: none)

    zunmbr('Q', 'L', 'C', M, NRHS, N, A, LDA, WORK(ITAUQ), B, LDB, WORK(IWORK),
        LWORK - IWORK + 1, INFO);

    // Generate right bidiagonalizing vectors in A
    // (CWorkspace: need 3*M, prefer 2*M+M*NB)
    // (RWorkspace: none)

    zungbr('P', M, N, M, A, LDA, WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1,
        INFO);
    IRWORK = IE + M;

    // Perform bidiagonal QR iteration,
    //    computing right singular vectors of A in A and
    //    multiplying B by transpose of left singular vectors
    // (CWorkspace: none)
    // (RWorkspace: need BDSPAC)

    zbdsqr('L', M, N, 0, NRHS, S, RWORK(IE), A, LDA, DUM.asMatrix(1), 1, B, LDB,
        RWORK(IRWORK), INFO);
    if (INFO.value != 0) {
      WORK[1] = MAXWRK.toComplex();
      return;
    }

    // Multiply B by reciprocals of singular values

    THR = max(RCOND * S[1], SFMIN);
    if (RCOND < ZERO) THR = max(EPS * S[1], SFMIN);
    RANK.value = 0;
    for (I = 1; I <= M; I++) {
      if (S[I] > THR) {
        zdrscl(NRHS, S[I], B(I, 1).asArray(), LDB);
        RANK.value++;
      } else {
        zlaset('F', 1, NRHS, Complex.zero, Complex.zero, B(I, 1), LDB);
      }
    }

    // Multiply B by right singular vectors of A
    // (CWorkspace: need N, prefer N*NRHS)
    // (RWorkspace: none)

    if (LWORK >= LDB * NRHS && NRHS > 1) {
      zgemm('C', 'N', N, NRHS, M, Complex.one, A, LDA, B, LDB, Complex.zero,
          WORK.asMatrix(LDB), LDB);
      zlacpy('G', N, NRHS, WORK.asMatrix(LDB), LDB, B, LDB);
    } else if (NRHS > 1) {
      CHUNK = LWORK ~/ N;
      for (I = 1; I <= NRHS; I += CHUNK) {
        BL = min(NRHS - I + 1, CHUNK);
        zgemm('C', 'N', N, BL, M, Complex.one, A, LDA, B(1, I), LDB,
            Complex.zero, WORK.asMatrix(N), N);
        zlacpy('F', N, BL, WORK.asMatrix(N), N, B(1, I), LDB);
      }
    } else if (NRHS == 1) {
      zgemv('C', M, N, Complex.one, A, LDA, B.asArray(), 1, Complex.zero, WORK,
          1);
      zcopy(N, WORK, 1, B.asArray(), 1);
    }
  }

  // Undo scaling

  if (IASCL == 1) {
    zlascl('G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO);
    dlascl('G', 0, 0, SMLNUM, ANRM, MINMN, 1, S.asMatrix(MINMN), MINMN, INFO);
  } else if (IASCL == 2) {
    zlascl('G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO);
    dlascl('G', 0, 0, BIGNUM, ANRM, MINMN, 1, S.asMatrix(MINMN), MINMN, INFO);
  }
  if (IBSCL == 1) {
    zlascl('G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO);
  } else if (IBSCL == 2) {
    zlascl('G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO);
  }
  WORK[1] = MAXWRK.toComplex();
}
