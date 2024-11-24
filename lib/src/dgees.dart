// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgebak.dart';
import 'package:lapack/src/dgebal.dart';
import 'package:lapack/src/dgehrd.dart';
import 'package:lapack/src/dhseqr.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dorghr.dart';
import 'package:lapack/src/dtrsen.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgees(
  final String JOBVS,
  final String SORT,
  final bool Function(double, double) SELECT,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Box<int> SDIM,
  final Array<double> WR_,
  final Array<double> WI_,
  final Matrix<double> VS_,
  final int LDVS,
  final Array<double> WORK_,
  final int LWORK,
  final Array<bool> BWORK_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final WR = WR_.having();
  final WI = WI_.having();
  final VS = VS_.having(ld: LDVS);
  final WORK = WORK_.having();
  final BWORK = BWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool CURSL, LASTSL, LQUERY, LST2SL, SCALEA, WANTST, WANTVS;
  int HSWORK, I, I1, I2, IBAL, INXT, IP, ITAU, IWRK, MAXWRK = 0, MINWRK;
  double ANRM, BIGNUM, CSCALE = 0, EPS, SMLNUM;
  final IDUM = Array<int>(1);
  final DUM = Array<double>(1);
  final IERR = Box(0),
      IEVAL = Box(0),
      ICOND = Box(0),
      IHI = Box(0),
      ILO = Box(0);
  final S = Box(0.0), SEP = Box(0.0);

  // Test the input arguments

  INFO.value = 0;
  LQUERY = (LWORK == -1);
  WANTVS = lsame(JOBVS, 'V');
  WANTST = lsame(SORT, 'S');
  if (!WANTVS && !lsame(JOBVS, 'N')) {
    INFO.value = -1;
  } else if (!WANTST && !lsame(SORT, 'N')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (LDA < max(1, N)) {
    INFO.value = -6;
  } else if (LDVS < 1 || (WANTVS && LDVS < N)) {
    INFO.value = -11;
  }

  // Compute workspace
  //  (Note: Comments in the code beginning "Workspace:" describe the
  //   minimal amount of workspace needed at that point in the code,
  //   as well as the preferred amount for good performance.
  //   NB refers to the optimal block size for the immediately
  //   following subroutine, as returned by ILAENV.
  //   HSWORK refers to the workspace preferred by DHSEQR, as
  //   calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
  //   the worst case.)

  if (INFO.value == 0) {
    if (N == 0) {
      MINWRK = 1;
      MAXWRK = 1;
    } else {
      MAXWRK = 2 * N + N * ilaenv(1, 'DGEHRD', ' ', N, 1, N, 0);
      MINWRK = 3 * N;

      dhseqr('S', JOBVS, N, 1, N, A, LDA, WR, WI, VS, LDVS, WORK, -1, IEVAL);
      HSWORK = WORK[1].toInt();

      if (!WANTVS) {
        MAXWRK = max(MAXWRK, N + HSWORK);
      } else {
        MAXWRK = max(
            MAXWRK, 2 * N + (N - 1) * ilaenv(1, 'DORGHR', ' ', N, 1, N, -1));
        MAXWRK = max(MAXWRK, N + HSWORK);
      }
    }
    WORK[1] = MAXWRK.toDouble();

    if (LWORK < MINWRK && !LQUERY) {
      INFO.value = -13;
    }
  }

  if (INFO.value != 0) {
    xerbla('DGEES', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) {
    SDIM.value = 0;
    return;
  }

  // Get machine constants

  EPS = dlamch('P');
  SMLNUM = dlamch('S');
  BIGNUM = ONE / SMLNUM;
  SMLNUM = sqrt(SMLNUM) / EPS;
  BIGNUM = ONE / SMLNUM;

  // Scale A if max element outside range [SMLNUM,BIGNUM]

  ANRM = dlange('M', N, N, A, LDA, DUM);
  SCALEA = false;
  if (ANRM > ZERO && ANRM < SMLNUM) {
    SCALEA = true;
    CSCALE = SMLNUM;
  } else if (ANRM > BIGNUM) {
    SCALEA = true;
    CSCALE = BIGNUM;
  }
  if (SCALEA) dlascl('G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR);

  // Permute the matrix to make it more nearly triangular
  // (Workspace: need N)

  IBAL = 1;
  dgebal('P', N, A, LDA, ILO, IHI, WORK(IBAL), IERR);

  // Reduce to upper Hessenberg form
  // (Workspace: need 3*N, prefer 2*N+N*NB)

  ITAU = N + IBAL;
  IWRK = N + ITAU;
  dgehrd(N, ILO.value, IHI.value, A, LDA, WORK(ITAU), WORK(IWRK),
      LWORK - IWRK + 1, IERR);

  if (WANTVS) {
    // Copy Householder vectors to VS

    dlacpy('L', N, N, A, LDA, VS, LDVS);

    // Generate orthogonal matrix in VS
    // (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)

    dorghr(N, ILO.value, IHI.value, VS, LDVS, WORK(ITAU), WORK(IWRK),
        LWORK - IWRK + 1, IERR);
  }

  SDIM.value = 0;

  // Perform QR iteration, accumulating Schur vectors in VS if desired
  // (Workspace: need N+1, prefer N+HSWORK (see comments) )

  IWRK = ITAU;
  dhseqr('S', JOBVS, N, ILO.value, IHI.value, A, LDA, WR, WI, VS, LDVS,
      WORK(IWRK), LWORK - IWRK + 1, IEVAL);
  if (IEVAL.value > 0) INFO.value = IEVAL.value;

  // Sort eigenvalues if desired

  if (WANTST && INFO.value == 0) {
    if (SCALEA) {
      dlascl('G', 0, 0, CSCALE, ANRM, N, 1, WR.asMatrix(N), N, IERR);
      dlascl('G', 0, 0, CSCALE, ANRM, N, 1, WI.asMatrix(N), N, IERR);
    }
    for (I = 1; I <= N; I++) {
      BWORK[I] = SELECT(WR[I], WI[I]);
    }

    // Reorder eigenvalues and transform Schur vectors
    // (Workspace: none needed)

    dtrsen('N', JOBVS, BWORK, N, A, LDA, VS, LDVS, WR, WI, SDIM, S, SEP,
        WORK(IWRK), LWORK - IWRK + 1, IDUM, 1, ICOND);
    if (ICOND.value > 0) INFO.value = N + ICOND.value;
  }

  if (WANTVS) {
    // Undo balancing
    // (Workspace: need N)

    dgebak('P', 'R', N, ILO.value, IHI.value, WORK(IBAL), N, VS, LDVS, IERR);
  }

  if (SCALEA) {
    // Undo scaling for the Schur form of A

    dlascl('H', 0, 0, CSCALE, ANRM, N, N, A, LDA, IERR);
    dcopy(N, A.asArray(), LDA + 1, WR, 1);
    if (CSCALE == SMLNUM) {
      // If scaling back towards underflow, adjust WI if an
      // offdiagonal element of a 2-by-2 block in the Schur form
      // underflows.

      if (IEVAL.value > 0) {
        I1 = IEVAL.value + 1;
        I2 = IHI.value - 1;
        final LWI = max(ILO.value - 1, 1);
        dlascl('G', 0, 0, CSCALE, ANRM, ILO.value - 1, 1, WI.asMatrix(LWI), LWI,
            IERR);
      } else if (WANTST) {
        I1 = 1;
        I2 = N - 1;
      } else {
        I1 = ILO.value;
        I2 = IHI.value - 1;
      }
      INXT = I1 - 1;
      for (I = I1; I <= I2; I++) {
        if (I < INXT) continue;
        if (WI[I] == ZERO) {
          INXT = I + 1;
        } else {
          if (A[I + 1][I] == ZERO) {
            WI[I] = ZERO;
            WI[I + 1] = ZERO;
          } else if (A[I + 1][I] != ZERO && A[I][I + 1] == ZERO) {
            WI[I] = ZERO;
            WI[I + 1] = ZERO;
            if (I > 1)
              // ignore: curly_braces_in_flow_control_structures
              dswap(I - 1, A(1, I).asArray(), 1, A(1, I + 1).asArray(), 1);
            if (N > I + 1) {
              dswap(N - I - 1, A(I, I + 2).asArray(), LDA,
                  A(I + 1, I + 2).asArray(), LDA);
            }
            if (WANTVS) {
              dswap(N, VS(1, I).asArray(), 1, VS(1, I + 1).asArray(), 1);
            }
            A[I][I + 1] = A[I + 1][I];
            A[I + 1][I] = ZERO;
          }
          INXT = I + 2;
        }
      }
    }

    // Undo scaling for the imaginary part of the eigenvalues
    final LWI = max(N - IEVAL.value, 1);
    dlascl('G', 0, 0, CSCALE, ANRM, N - IEVAL.value, 1,
        WI(IEVAL.value + 1).asMatrix(LWI), LWI, IERR);
  }

  if (WANTST && INFO.value == 0) {
    // Check if reordering successful

    LASTSL = true;
    LST2SL = true;
    SDIM.value = 0;
    IP = 0;
    for (I = 1; I <= N; I++) {
      CURSL = SELECT(WR[I], WI[I]);
      if (WI[I] == ZERO) {
        if (CURSL) SDIM.value++;
        IP = 0;
        if (CURSL && !LASTSL) INFO.value = N + 2;
      } else {
        if (IP == 1) {
          // Last eigenvalue of conjugate pair

          CURSL = CURSL || LASTSL;
          LASTSL = CURSL;
          if (CURSL) SDIM.value += 2;
          IP = -1;
          if (CURSL && !LST2SL) INFO.value = N + 2;
        } else {
          // First eigenvalue of conjugate pair

          IP = 1;
        }
      }
      LST2SL = LASTSL;
      LASTSL = CURSL;
    }
  }

  WORK[1] = MAXWRK.toDouble();
}
