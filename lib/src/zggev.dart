// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zgeqrf.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zggbak.dart';
import 'package:dart_lapack/src/zggbal.dart';
import 'package:dart_lapack/src/zgghrd.dart';
import 'package:dart_lapack/src/zhgeqz.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlascl.dart';
import 'package:dart_lapack/src/zlaset.dart';
import 'package:dart_lapack/src/ztgevc.dart';
import 'package:dart_lapack/src/zungqr.dart';
import 'package:dart_lapack/src/zunmqr.dart';

void zggev(
  final String JOBVL,
  final String JOBVR,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<Complex> ALPHA_,
  final Array<Complex> BETA_,
  final Matrix<Complex> VL_,
  final int LDVL,
  final Matrix<Complex> VR_,
  final int LDVR,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final VL = VL_.having(ld: LDVL);
  final VR = VR_.having(ld: LDVR);
  final ALPHA = ALPHA_.having();
  final BETA = BETA_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool ILASCL, ILBSCL, ILV, ILVL, ILVR, LQUERY;
  String CHTEMP;
  int ICOLS,
      IJOBVL,
      IJOBVR,
      ILEFT,
      IRIGHT,
      IROWS,
      IRWRK,
      ITAU,
      IWRK,
      JC,
      JR,
      LWKMIN,
      LWKOPT = 0;
  double ANRM, ANRMTO = 0, BIGNUM, BNRM, BNRMTO = 0, EPS, SMLNUM, TEMP;
  final LDUMMA = Array<bool>(1);
  final IERR = Box(0), IHI = Box(0), ILO = Box(0), IN = Box(0);

  // Decode the input arguments

  if (lsame(JOBVL, 'N')) {
    IJOBVL = 1;
    ILVL = false;
  } else if (lsame(JOBVL, 'V')) {
    IJOBVL = 2;
    ILVL = true;
  } else {
    IJOBVL = -1;
    ILVL = false;
  }

  if (lsame(JOBVR, 'N')) {
    IJOBVR = 1;
    ILVR = false;
  } else if (lsame(JOBVR, 'V')) {
    IJOBVR = 2;
    ILVR = true;
  } else {
    IJOBVR = -1;
    ILVR = false;
  }
  ILV = ILVL || ILVR;

  // Test the input arguments

  INFO.value = 0;
  LQUERY = (LWORK == -1);
  if (IJOBVL <= 0) {
    INFO.value = -1;
  } else if (IJOBVR <= 0) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDB < max(1, N)) {
    INFO.value = -7;
  } else if (LDVL < 1 || (ILVL && LDVL < N)) {
    INFO.value = -11;
  } else if (LDVR < 1 || (ILVR && LDVR < N)) {
    INFO.value = -13;
  }

  // Compute workspace
  //  (Note: Comments in the code beginning "Workspace:" describe the
  //   minimal amount of workspace needed at that point in the code,
  //   as well as the preferred amount for good performance.
  //   NB refers to the optimal block size for the immediately
  //   following subroutine, as returned by ILAENV. The workspace is
  //   computed assuming ILO = 1 and IHI = N, the worst case.)

  if (INFO.value == 0) {
    LWKMIN = max(1, 2 * N);
    LWKOPT = max(1, N + N * ilaenv(1, 'ZGEQRF', ' ', N, 1, N, 0));
    LWKOPT = max(LWKOPT, N + N * ilaenv(1, 'ZUNMQR', ' ', N, 1, N, 0));
    if (ILVL) {
      LWKOPT = max(LWKOPT, N + N * ilaenv(1, 'ZUNGQR', ' ', N, 1, N, -1));
    }
    WORK[1] = LWKOPT.toComplex();

    if (LWORK < LWKMIN && !LQUERY) INFO.value = -15;
  }

  if (INFO.value != 0) {
    xerbla('ZGGEV', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Get machine constants

  EPS = dlamch('E') * dlamch('B');
  SMLNUM = dlamch('S');
  BIGNUM = ONE / SMLNUM;
  SMLNUM = sqrt(SMLNUM) / EPS;
  BIGNUM = ONE / SMLNUM;

  // Scale A if max element outside range [SMLNUM,BIGNUM]

  ANRM = zlange('M', N, N, A, LDA, RWORK);
  ILASCL = false;
  if (ANRM > ZERO && ANRM < SMLNUM) {
    ANRMTO = SMLNUM;
    ILASCL = true;
  } else if (ANRM > BIGNUM) {
    ANRMTO = BIGNUM;
    ILASCL = true;
  }
  if (ILASCL) zlascl('G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR);

  // Scale B if max element outside range [SMLNUM,BIGNUM]

  BNRM = zlange('M', N, N, B, LDB, RWORK);
  ILBSCL = false;
  if (BNRM > ZERO && BNRM < SMLNUM) {
    BNRMTO = SMLNUM;
    ILBSCL = true;
  } else if (BNRM > BIGNUM) {
    BNRMTO = BIGNUM;
    ILBSCL = true;
  }
  if (ILBSCL) zlascl('G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR);

  // Permute the matrices A, B to isolate eigenvalues if possible
  // (Real Workspace: need 6*N)

  ILEFT = 1;
  IRIGHT = N + 1;
  IRWRK = IRIGHT + N;
  zggbal('P', N, A, LDA, B, LDB, ILO, IHI, RWORK(ILEFT), RWORK(IRIGHT),
      RWORK(IRWRK), IERR);

  // Reduce B to triangular form (QR decomposition of B)
  // (Complex Workspace: need N, prefer N*NB)

  IROWS = IHI.value + 1 - ILO.value;
  if (ILV) {
    ICOLS = N + 1 - ILO.value;
  } else {
    ICOLS = IROWS;
  }
  ITAU = 1;
  IWRK = ITAU + IROWS;
  zgeqrf(IROWS, ICOLS, B(ILO.value, ILO.value), LDB, WORK(ITAU), WORK(IWRK),
      LWORK + 1 - IWRK, IERR);

  // Apply the orthogonal transformation to matrix A
  // (Complex Workspace: need N, prefer N*NB)

  zunmqr(
      'L',
      'C',
      IROWS,
      ICOLS,
      IROWS,
      B(ILO.value, ILO.value),
      LDB,
      WORK(ITAU),
      A(ILO.value, ILO.value),
      LDA,
      WORK(IWRK),
      LWORK + 1 - IWRK,
      IERR);

  // Initialize VL
  // (Complex Workspace: need N, prefer N*NB)

  if (ILVL) {
    zlaset('Full', N, N, Complex.zero, Complex.one, VL, LDVL);
    if (IROWS > 1) {
      zlacpy('L', IROWS - 1, IROWS - 1, B(ILO.value + 1, ILO.value), LDB,
          VL(ILO.value + 1, ILO.value), LDVL);
    }
    zungqr(IROWS, IROWS, IROWS, VL(ILO.value, ILO.value), LDVL, WORK(ITAU),
        WORK(IWRK), LWORK + 1 - IWRK, IERR);
  }

  // Initialize VR

  if (ILVR) zlaset('Full', N, N, Complex.zero, Complex.one, VR, LDVR);

  // Reduce to generalized Hessenberg form

  if (ILV) {
    // Eigenvectors requested -- work on whole matrix.

    zgghrd(JOBVL, JOBVR, N, ILO.value, IHI.value, A, LDA, B, LDB, VL, LDVL, VR,
        LDVR, IERR);
  } else {
    zgghrd('N', 'N', IROWS, 1, IROWS, A(ILO.value, ILO.value), LDA,
        B(ILO.value, ILO.value), LDB, VL, LDVL, VR, LDVR, IERR);
  }

  // Perform QZ algorithm (Compute eigenvalues, and optionally, the
  // Schur form and Schur vectors)
  // (Complex Workspace: need N)
  // (Real Workspace: need N)

  IWRK = ITAU;
  if (ILV) {
    CHTEMP = 'S';
  } else {
    CHTEMP = 'E';
  }
  zhgeqz(
      CHTEMP,
      JOBVL,
      JOBVR,
      N,
      ILO.value,
      IHI.value,
      A,
      LDA,
      B,
      LDB,
      ALPHA,
      BETA,
      VL,
      LDVL,
      VR,
      LDVR,
      WORK(IWRK),
      LWORK + 1 - IWRK,
      RWORK(IRWRK),
      IERR);
  if (IERR.value != 0) {
    if (IERR.value > 0 && IERR.value <= N) {
      INFO.value = IERR.value;
    } else if (IERR.value > N && IERR.value <= 2 * N) {
      INFO.value = IERR.value - N;
    } else {
      INFO.value = N + 1;
    }
  } else

  // Compute Eigenvectors
  // (Real Workspace: need 2*N)
  // (Complex Workspace: need 2*N)

  if (ILV) {
    if (ILVL) {
      if (ILVR) {
        CHTEMP = 'B';
      } else {
        CHTEMP = 'L';
      }
    } else {
      CHTEMP = 'R';
    }

    ztgevc(CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, N, IN,
        WORK(IWRK), RWORK(IRWRK), IERR);
    if (IERR.value != 0) {
      INFO.value = N + 2;
    } else {
      // Undo balancing on VL and VR and normalization
      // (Workspace: none needed)

      if (ILVL) {
        zggbak('P', 'L', N, ILO.value, IHI.value, RWORK(ILEFT), RWORK(IRIGHT),
            N, VL, LDVL, IERR);
        for (JC = 1; JC <= N; JC++) {
          TEMP = ZERO;
          for (JR = 1; JR <= N; JR++) {
            TEMP = max(TEMP, VL[JR][JC].cabs1());
          }
          if (TEMP < SMLNUM) continue;
          TEMP = ONE / TEMP;
          for (JR = 1; JR <= N; JR++) {
            VL[JR][JC] *= TEMP.toComplex();
          }
        }
      }
      if (ILVR) {
        zggbak('P', 'R', N, ILO.value, IHI.value, RWORK(ILEFT), RWORK(IRIGHT),
            N, VR, LDVR, IERR);
        for (JC = 1; JC <= N; JC++) {
          TEMP = ZERO;
          for (JR = 1; JR <= N; JR++) {
            TEMP = max(TEMP, VR[JR][JC].cabs1());
          }
          if (TEMP < SMLNUM) continue;
          TEMP = ONE / TEMP;
          for (JR = 1; JR <= N; JR++) {
            VR[JR][JC] *= TEMP.toComplex();
          }
        }
      }
    }
  }

  // Undo scaling if necessary

  if (ILASCL) zlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA.asMatrix(N), N, IERR);

  if (ILBSCL) zlascl('G', 0, 0, BNRMTO, BNRM, N, 1, BETA.asMatrix(N), N, IERR);

  WORK[1] = LWKOPT.toComplex();
}
