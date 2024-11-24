// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/blas/drot.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgebak.dart';
import 'package:lapack/src/dgebal.dart';
import 'package:lapack/src/dgehrd.dart';
import 'package:lapack/src/dhseqr.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlapy2.dart';
import 'package:lapack/src/dlartg.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dorghr.dart';
import 'package:lapack/src/dtrevc3.dart';
import 'package:lapack/src/dtrsna.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgeevx(
  final String BALANC,
  final String JOBVL,
  final String JOBVR,
  final String SENSE,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> WR_,
  final Array<double> WI_,
  final Matrix<double> VL_,
  final int LDVL,
  final Matrix<double> VR_,
  final int LDVR,
  final Box<int> ILO,
  final Box<int> IHI,
  final Array<double> SCALE_,
  final Box<double> ABNRM,
  final Array<double> RCONDE_,
  final Array<double> RCONDV_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final WR = WR_.having();
  final WI = WI_.having();
  final VL = VL_.having(ld: LDVL);
  final VR = VR_.having(ld: LDVR);
  final SCALE = SCALE_.having();
  final RCONDE = RCONDE_.having();
  final RCONDV = RCONDV_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool LQUERY, SCALEA, WANTVL, WANTVR, WNTSNB, WNTSNE, WNTSNN, WNTSNV;
  String JOB, SIDE = '';
  int HSWORK, I, ITAU, IWRK, K, LWORK_TREVC, MAXWRK = 0, MINWRK;
  double ANRM, BIGNUM, CSCALE = 0, EPS, SCL, SMLNUM;
  final SELECT = Array<bool>(1);
  final DUM = Array<double>(1);
  final IERR = Box(0), NOUT = Box(0), ICOND = Box(0);
  final CS = Box(0.0), SN = Box(0.0), R = Box(0.0);

  // Test the input arguments

  INFO.value = 0;
  LQUERY = (LWORK == -1);
  WANTVL = lsame(JOBVL, 'V');
  WANTVR = lsame(JOBVR, 'V');
  WNTSNN = lsame(SENSE, 'N');
  WNTSNE = lsame(SENSE, 'E');
  WNTSNV = lsame(SENSE, 'V');
  WNTSNB = lsame(SENSE, 'B');
  if (!(lsame(BALANC, 'N') ||
      lsame(BALANC, 'S') ||
      lsame(BALANC, 'P') ||
      lsame(BALANC, 'B'))) {
    INFO.value = -1;
  } else if (!WANTVL && !lsame(JOBVL, 'N')) {
    INFO.value = -2;
  } else if (!WANTVR && !lsame(JOBVR, 'N')) {
    INFO.value = -3;
  } else if (!(WNTSNN || WNTSNE || WNTSNB || WNTSNV) ||
      ((WNTSNE || WNTSNB) && !(WANTVL && WANTVR))) {
    INFO.value = -4;
  } else if (N < 0) {
    INFO.value = -5;
  } else if (LDA < max(1, N)) {
    INFO.value = -7;
  } else if (LDVL < 1 || (WANTVL && LDVL < N)) {
    INFO.value = -11;
  } else if (LDVR < 1 || (WANTVR && LDVR < N)) {
    INFO.value = -13;
  }

  // Compute workspace
  // (Note: Comments in the code beginning "Workspace:" describe the
  // minimal amount of workspace needed at that point in the code,
  // as well as the preferred amount for good performance.
  // NB refers to the optimal block size for the immediately
  // following subroutine, as returned by ILAENV.
  // HSWORK refers to the workspace preferred by DHSEQR, as
  // calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
  // the worst case.)

  if (INFO.value == 0) {
    if (N == 0) {
      MINWRK = 1;
      MAXWRK = 1;
    } else {
      MAXWRK = N + N * ilaenv(1, 'DGEHRD', ' ', N, 1, N, 0);

      if (WANTVL) {
        dtrevc3('L', 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK,
            -1, IERR);
        LWORK_TREVC = WORK[1].toInt();
        MAXWRK = max(MAXWRK, N + LWORK_TREVC);
        dhseqr('S', 'V', N, 1, N, A, LDA, WR, WI, VL, LDVL, WORK, -1, INFO);
      } else if (WANTVR) {
        dtrevc3('R', 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK,
            -1, IERR);
        LWORK_TREVC = WORK[1].toInt();
        MAXWRK = max(MAXWRK, N + LWORK_TREVC);
        dhseqr('S', 'V', N, 1, N, A, LDA, WR, WI, VR, LDVR, WORK, -1, INFO);
      } else {
        if (WNTSNN) {
          dhseqr('E', 'N', N, 1, N, A, LDA, WR, WI, VR, LDVR, WORK, -1, INFO);
        } else {
          dhseqr('S', 'N', N, 1, N, A, LDA, WR, WI, VR, LDVR, WORK, -1, INFO);
        }
      }
      HSWORK = WORK[1].toInt();

      if (!WANTVL && !WANTVR) {
        MINWRK = 2 * N;
        if (!WNTSNN) MINWRK = max(MINWRK, N * N + 6 * N);
        MAXWRK = max(MAXWRK, HSWORK);
        if (!WNTSNN) MAXWRK = max(MAXWRK, N * N + 6 * N);
      } else {
        MINWRK = 3 * N;
        if (!WNTSNN && !WNTSNE) MINWRK = max(MINWRK, N * N + 6 * N);
        MAXWRK = max(MAXWRK, HSWORK);
        MAXWRK =
            max(MAXWRK, N + (N - 1) * ilaenv(1, 'DORGHR', ' ', N, 1, N, -1));
        if (!WNTSNN && !WNTSNE) MAXWRK = max(MAXWRK, N * N + 6 * N);
        MAXWRK = max(MAXWRK, 3 * N);
      }
      MAXWRK = max(MAXWRK, MINWRK);
    }
    WORK[1] = MAXWRK.toDouble();

    if (LWORK < MINWRK && !LQUERY) {
      INFO.value = -21;
    }
  }

  if (INFO.value != 0) {
    xerbla('DGEEVX', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Get machine constants

  EPS = dlamch('P');
  SMLNUM = dlamch('S');
  BIGNUM = ONE / SMLNUM;
  SMLNUM = sqrt(SMLNUM) / EPS;
  BIGNUM = ONE / SMLNUM;

  // Scale A if max element outside range [SMLNUM,BIGNUM]

  ICOND.value = 0;
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

  // Balance the matrix and compute ABNRM

  dgebal(BALANC, N, A, LDA, ILO, IHI, SCALE, IERR);
  ABNRM.value = dlange('1', N, N, A, LDA, DUM);
  if (SCALEA) {
    DUM[1] = ABNRM.value;
    dlascl('G', 0, 0, CSCALE, ANRM, 1, 1, DUM.asMatrix(1), 1, IERR);
    ABNRM.value = DUM[1];
  }

  // Reduce to upper Hessenberg form
  // (Workspace: need 2*N, prefer N+N*NB)

  ITAU = 1;
  IWRK = ITAU + N;
  dgehrd(N, ILO.value, IHI.value, A, LDA, WORK(ITAU), WORK(IWRK),
      LWORK - IWRK + 1, IERR);

  if (WANTVL) {
    // Want left eigenvectors
    // Copy Householder vectors to VL

    SIDE = 'L';
    dlacpy('L', N, N, A, LDA, VL, LDVL);

    // Generate orthogonal matrix in VL
    // (Workspace: need 2*N-1, prefer N+(N-1)*NB)

    dorghr(N, ILO.value, IHI.value, VL, LDVL, WORK(ITAU), WORK(IWRK),
        LWORK - IWRK + 1, IERR);

    // Perform QR iteration, accumulating Schur vectors in VL
    // (Workspace: need 1, prefer HSWORK (see comments) )

    IWRK = ITAU;
    dhseqr('S', 'V', N, ILO.value, IHI.value, A, LDA, WR, WI, VL, LDVL,
        WORK(IWRK), LWORK - IWRK + 1, INFO);

    if (WANTVR) {
      // Want left and right eigenvectors
      // Copy Schur vectors to VR

      SIDE = 'B';
      dlacpy('F', N, N, VL, LDVL, VR, LDVR);
    }
  } else if (WANTVR) {
    // Want right eigenvectors
    // Copy Householder vectors to VR

    SIDE = 'R';
    dlacpy('L', N, N, A, LDA, VR, LDVR);

    // Generate orthogonal matrix in VR
    // (Workspace: need 2*N-1, prefer N+(N-1)*NB)

    dorghr(N, ILO.value, IHI.value, VR, LDVR, WORK(ITAU), WORK(IWRK),
        LWORK - IWRK + 1, IERR);

    // Perform QR iteration, accumulating Schur vectors in VR
    // (Workspace: need 1, prefer HSWORK (see comments) )

    IWRK = ITAU;
    dhseqr('S', 'V', N, ILO.value, IHI.value, A, LDA, WR, WI, VR, LDVR,
        WORK(IWRK), LWORK - IWRK + 1, INFO);
  } else {
    // Compute eigenvalues only
    // If condition numbers desired, compute Schur form

    if (WNTSNN) {
      JOB = 'E';
    } else {
      JOB = 'S';
    }

    // (Workspace: need 1, prefer HSWORK (see comments) )

    IWRK = ITAU;
    dhseqr(JOB, 'N', N, ILO.value, IHI.value, A, LDA, WR, WI, VR, LDVR,
        WORK(IWRK), LWORK - IWRK + 1, INFO);
  }

  // If INFO != 0 from DHSEQR, then quit

  if (INFO.value == 0) {
    if (WANTVL || WANTVR) {
      // Compute left and/or right eigenvectors
      // (Workspace: need 3*N, prefer N + 2*N*NB)

      dtrevc3(SIDE, 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT,
          WORK(IWRK), LWORK - IWRK + 1, IERR);
    }

    // Compute condition numbers if desired
    // (Workspace: need N*N+6*N unless SENSE = 'E')

    if (!WNTSNN) {
      dtrsna(SENSE, 'A', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, RCONDE, RCONDV,
          N, NOUT, WORK(IWRK).asMatrix(N), N, IWORK, ICOND);
    }

    if (WANTVL) {
      // Undo balancing of left eigenvectors

      dgebak(BALANC, 'L', N, ILO.value, IHI.value, SCALE, N, VL, LDVL, IERR);

      // Normalize left eigenvectors and make largest component real

      for (I = 1; I <= N; I++) {
        if (WI[I] == ZERO) {
          SCL = ONE / dnrm2(N, VL(1, I).asArray(), 1);
          dscal(N, SCL, VL(1, I).asArray(), 1);
        } else if (WI[I] > ZERO) {
          SCL = ONE /
              dlapy2(dnrm2(N, VL(1, I).asArray(), 1),
                  dnrm2(N, VL(1, I + 1).asArray(), 1));
          dscal(N, SCL, VL(1, I).asArray(), 1);
          dscal(N, SCL, VL(1, I + 1).asArray(), 1);
          for (K = 1; K <= N; K++) {
            WORK[K] = pow(VL[K][I], 2).toDouble() + pow(VL[K][I + 1], 2);
          }
          K = idamax(N, WORK, 1);
          dlartg(VL[K][I], VL[K][I + 1], CS, SN, R);
          drot(N, VL(1, I).asArray(), 1, VL(1, I + 1).asArray(), 1, CS.value,
              SN.value);
          VL[K][I + 1] = ZERO;
        }
      }
    }

    if (WANTVR) {
      // Undo balancing of right eigenvectors

      dgebak(BALANC, 'R', N, ILO.value, IHI.value, SCALE, N, VR, LDVR, IERR);

      // Normalize right eigenvectors and make largest component real

      for (I = 1; I <= N; I++) {
        if (WI[I] == ZERO) {
          SCL = ONE / dnrm2(N, VR(1, I).asArray(), 1);
          dscal(N, SCL, VR(1, I).asArray(), 1);
        } else if (WI[I] > ZERO) {
          SCL = ONE /
              dlapy2(dnrm2(N, VR(1, I).asArray(), 1),
                  dnrm2(N, VR(1, I + 1).asArray(), 1));
          dscal(N, SCL, VR(1, I).asArray(), 1);
          dscal(N, SCL, VR(1, I + 1).asArray(), 1);
          for (K = 1; K <= N; K++) {
            WORK[K] = pow(VR[K][I], 2).toDouble() + pow(VR[K][I + 1], 2);
          }
          K = idamax(N, WORK, 1);
          dlartg(VR[K][I], VR[K][I + 1], CS, SN, R);
          drot(N, VR(1, I).asArray(), 1, VR(1, I + 1).asArray(), 1, CS.value,
              SN.value);
          VR[K][I + 1] = ZERO;
        }
      }
    }
  }

  // Undo scaling if necessary

  if (SCALEA) {
    final ld = max(N - INFO.value, 1);
    dlascl('G', 0, 0, CSCALE, ANRM, N - INFO.value, 1,
        WR(INFO.value + 1).asMatrix(ld), ld, IERR);
    dlascl('G', 0, 0, CSCALE, ANRM, N - INFO.value, 1,
        WI(INFO.value + 1).asMatrix(ld), ld, IERR);
    if (INFO.value == 0) {
      if ((WNTSNV || WNTSNB) && ICOND.value == 0) {
        dlascl('G', 0, 0, CSCALE, ANRM, N, 1, RCONDV.asMatrix(N), N, IERR);
      }
    } else {
      dlascl(
          'G', 0, 0, CSCALE, ANRM, ILO.value - 1, 1, WR.asMatrix(N), N, IERR);
      dlascl(
          'G', 0, 0, CSCALE, ANRM, ILO.value - 1, 1, WI.asMatrix(N), N, IERR);
    }
  }

  WORK[1] = MAXWRK.toDouble();
}
