// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgeqrf.dart';
import 'package:dart_lapack/src/dggbak.dart';
import 'package:dart_lapack/src/dggbal.dart';
import 'package:dart_lapack/src/dgghrd.dart';
import 'package:dart_lapack/src/dhgeqz.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dlascl.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dorgqr.dart';
import 'package:dart_lapack/src/dormqr.dart';
import 'package:dart_lapack/src/dtgevc.dart';
import 'package:dart_lapack/src/dtgsna.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dggevx(
  final String BALANC,
  final String JOBVL,
  final String JOBVR,
  final String SENSE,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Array<double> ALPHAR_,
  final Array<double> ALPHAI_,
  final Array<double> BETA_,
  final Matrix<double> VL_,
  final int LDVL,
  final Matrix<double> VR_,
  final int LDVR,
  final Box<int> ILO,
  final Box<int> IHI,
  final Array<double> LSCALE_,
  final Array<double> RSCALE_,
  final Box<double> ABNRM,
  final Box<double> BBNRM,
  final Array<double> RCONDE_,
  final Array<double> RCONDV_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final Array<bool> BWORK_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final ALPHAR = ALPHAR_.having();
  final ALPHAI = ALPHAI_.having();
  final BETA = BETA_.having();
  final VL = VL_.having(ld: LDVL);
  final VR = VR_.having(ld: LDVR);
  final LSCALE = LSCALE_.having();
  final RSCALE = RSCALE_.having();
  final RCONDE = RCONDE_.having();
  final RCONDV = RCONDV_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  final BWORK = BWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool ILASCL,
      ILBSCL,
      ILV,
      ILVL,
      ILVR,
      LQUERY,
      NOSCL,
      PAIR,
      WANTSB,
      WANTSE,
      WANTSN,
      WANTSV;
  String CHTEMP;
  int I,
      ICOLS,
      IJOBVL,
      IJOBVR,
      IROWS,
      ITAU,
      IWRK,
      IWRK1,
      J,
      JC,
      JR,
      MAXWRK = 0,
      MINWRK,
      MM;
  double ANRM, ANRMTO = 0, BIGNUM, BNRM, BNRMTO = 0, EPS, SMLNUM, TEMP;
  final IERR = Box(0), IN = Box(0), M = Box(0);
  final LDUMMA = Array<bool>(1);

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

  NOSCL = lsame(BALANC, 'N') || lsame(BALANC, 'P');
  WANTSN = lsame(SENSE, 'N');
  WANTSE = lsame(SENSE, 'E');
  WANTSV = lsame(SENSE, 'V');
  WANTSB = lsame(SENSE, 'B');

  // Test the input arguments

  INFO.value = 0;
  LQUERY = (LWORK == -1);
  if (!(lsame(BALANC, 'N') ||
      lsame(BALANC, 'S') ||
      lsame(BALANC, 'P') ||
      lsame(BALANC, 'B'))) {
    INFO.value = -1;
  } else if (IJOBVL <= 0) {
    INFO.value = -2;
  } else if (IJOBVR <= 0) {
    INFO.value = -3;
  } else if (!(WANTSN || WANTSE || WANTSB || WANTSV)) {
    INFO.value = -4;
  } else if (N < 0) {
    INFO.value = -5;
  } else if (LDA < max(1, N)) {
    INFO.value = -7;
  } else if (LDB < max(1, N)) {
    INFO.value = -9;
  } else if (LDVL < 1 || (ILVL && LDVL < N)) {
    INFO.value = -14;
  } else if (LDVR < 1 || (ILVR && LDVR < N)) {
    INFO.value = -16;
  }

  // Compute workspace
  // (Note: Comments in the code beginning "Workspace:" describe the
  // minimal amount of workspace needed at that point in the code,
  // as well as the preferred amount for good performance.
  // NB refers to the optimal block size for the immediately
  // following subroutine, as returned by ilaenv. The workspace is
  // computed assuming ILO = 1 and IHI = N, the worst case.)

  if (INFO.value == 0) {
    if (N == 0) {
      MINWRK = 1;
      MAXWRK = 1;
    } else {
      if (NOSCL && !ILV) {
        MINWRK = 2 * N;
      } else {
        MINWRK = 6 * N;
      }
      if (WANTSE || WANTSB) {
        MINWRK = 10 * N;
      }
      if (WANTSV || WANTSB) {
        MINWRK = max(MINWRK, 2 * N * (N + 4) + 16);
      }
      MAXWRK = MINWRK;
      MAXWRK = max(MAXWRK, N + N * ilaenv(1, 'DGEQRF', ' ', N, 1, N, 0));
      MAXWRK = max(MAXWRK, N + N * ilaenv(1, 'DORMQR', ' ', N, 1, N, 0));
      if (ILVL) {
        MAXWRK = max(MAXWRK, N + N * ilaenv(1, 'DORGQR', ' ', N, 1, N, 0));
      }
    }
    WORK[1] = MAXWRK.toDouble();

    if (LWORK < MINWRK && !LQUERY) {
      INFO.value = -26;
    }
  }

  if (INFO.value != 0) {
    xerbla('DGGEVX', -INFO.value);
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

  ANRM = dlange('M', N, N, A, LDA, WORK);
  ILASCL = false;
  if (ANRM > ZERO && ANRM < SMLNUM) {
    ANRMTO = SMLNUM;
    ILASCL = true;
  } else if (ANRM > BIGNUM) {
    ANRMTO = BIGNUM;
    ILASCL = true;
  }
  if (ILASCL) dlascl('G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR);

  // Scale B if max element outside range [SMLNUM,BIGNUM]

  BNRM = dlange('M', N, N, B, LDB, WORK);
  ILBSCL = false;
  if (BNRM > ZERO && BNRM < SMLNUM) {
    BNRMTO = SMLNUM;
    ILBSCL = true;
  } else if (BNRM > BIGNUM) {
    BNRMTO = BIGNUM;
    ILBSCL = true;
  }
  if (ILBSCL) dlascl('G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR);

  // Permute and/or balance the matrix pair (A,B)
  // (Workspace: need 6*N if BALANC = 'S' or 'B', 1 otherwise)

  dggbal(BALANC, N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, WORK, IERR);

  // Compute ABNRM and BBNRM

  ABNRM.value = dlange('1', N, N, A, LDA, WORK(1));
  if (ILASCL) {
    WORK[1] = ABNRM.value;
    dlascl('G', 0, 0, ANRMTO, ANRM, 1, 1, WORK.asMatrix(1), 1, IERR);
    ABNRM.value = WORK[1];
  }

  BBNRM.value = dlange('1', N, N, B, LDB, WORK(1));
  if (ILBSCL) {
    WORK[1] = BBNRM.value;
    dlascl('G', 0, 0, BNRMTO, BNRM, 1, 1, WORK.asMatrix(1), 1, IERR);
    BBNRM.value = WORK[1];
  }

  // Reduce B to triangular form (QR decomposition of B)
  // (Workspace: need N, prefer N*NB )

  IROWS = IHI.value + 1 - ILO.value;
  if (ILV || !WANTSN) {
    ICOLS = N + 1 - ILO.value;
  } else {
    ICOLS = IROWS;
  }
  ITAU = 1;
  IWRK = ITAU + IROWS;
  dgeqrf(IROWS, ICOLS, B(ILO.value, ILO.value), LDB, WORK(ITAU), WORK(IWRK),
      LWORK + 1 - IWRK, IERR);

  // Apply the orthogonal transformation to A
  // (Workspace: need N, prefer N*NB)

  dormqr(
      'L',
      'T',
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

  // Initialize VL and/or VR
  // (Workspace: need N, prefer N*NB)

  if (ILVL) {
    dlaset('Full', N, N, ZERO, ONE, VL, LDVL);
    if (IROWS > 1) {
      dlacpy('L', IROWS - 1, IROWS - 1, B(ILO.value + 1, ILO.value), LDB,
          VL(ILO.value + 1, ILO.value), LDVL);
    }
    dorgqr(IROWS, IROWS, IROWS, VL(ILO.value, ILO.value), LDVL, WORK(ITAU),
        WORK(IWRK), LWORK + 1 - IWRK, IERR);
  }

  if (ILVR) dlaset('Full', N, N, ZERO, ONE, VR, LDVR);

  // Reduce to generalized Hessenberg form
  // (Workspace: none needed)

  if (ILV || !WANTSN) {
    // Eigenvectors requested -- work on whole matrix.

    dgghrd(JOBVL, JOBVR, N, ILO.value, IHI.value, A, LDA, B, LDB, VL, LDVL, VR,
        LDVR, IERR);
  } else {
    dgghrd('N', 'N', IROWS, 1, IROWS, A(ILO.value, ILO.value), LDA,
        B(ILO.value, ILO.value), LDB, VL, LDVL, VR, LDVR, IERR);
  }

  // Perform QZ algorithm (Compute eigenvalues, and optionally, the
  // Schur forms and Schur vectors)
  // (Workspace: need N)

  if (ILV || !WANTSN) {
    CHTEMP = 'S';
  } else {
    CHTEMP = 'E';
  }
  try {
    dhgeqz(CHTEMP, JOBVL, JOBVR, N, ILO.value, IHI.value, A, LDA, B, LDB,
        ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, IERR);
    if (IERR.value != 0) {
      if (IERR.value > 0 && IERR.value <= N) {
        INFO.value = IERR.value;
      } else if (IERR.value > N && IERR.value <= 2 * N) {
        INFO.value = IERR.value - N;
      } else {
        INFO.value = N + 1;
      }
      return;
    }

    // Compute Eigenvectors and estimate condition numbers if desired
    // (Workspace: DTGEVC: need 6*N
    //             DTGSNA: need 2*N*(N+2)+16 if SENSE = 'V' or 'B',
    //                     need N otherwise )
    if (ILV || !WANTSN) {
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

        dtgevc(CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, N,
            IN, WORK, IERR);
        if (IERR.value != 0) {
          INFO.value = N + 2;
          return;
        }
      }

      if (!WANTSN) {
        // compute eigenvectors (DTGEVC) and estimate condition
        // numbers (DTGSNA). Note that the definition of the condition
        // number is not invariant under transformation (u,v) to
        // (Q*u, Z*v), where (u,v) are eigenvectors of the generalized
        // Schur form (S,T), Q and Z are orthogonal matrices. In order
        // to avoid using extra 2*N*N workspace, we have to recalculate
        // eigenvectors and estimate one condition numbers at a time.

        PAIR = false;
        for (I = 1; I <= N; I++) {
          if (PAIR) {
            PAIR = false;
            continue;
          }
          MM = 1;
          if (I < N) {
            if (A[I + 1][I] != ZERO) {
              PAIR = true;
              MM = 2;
            }
          }

          for (J = 1; J <= N; J++) {
            BWORK[J] = false;
          }
          if (MM == 1) {
            BWORK[I] = true;
          } else if (MM == 2) {
            BWORK[I] = true;
            BWORK[I + 1] = true;
          }

          IWRK = MM * N + 1;
          IWRK1 = IWRK + MM * N;

          // Compute a pair of left and right eigenvectors.
          // (compute workspace: need up to 4*N + 6*N)

          if (WANTSE || WANTSB) {
            dtgevc('B', 'S', BWORK, N, A, LDA, B, LDB, WORK.asMatrix(N), N,
                WORK(IWRK).asMatrix(N), N, MM, M, WORK(IWRK1), IERR);
            if (IERR.value != 0) {
              INFO.value = N + 2;
              return;
            }
          }

          dtgsna(
              SENSE,
              'S',
              BWORK,
              N,
              A,
              LDA,
              B,
              LDB,
              WORK.asMatrix(N),
              N,
              WORK(IWRK).asMatrix(N),
              N,
              RCONDE(I),
              RCONDV(I),
              MM,
              M,
              WORK(IWRK1),
              LWORK - IWRK1 + 1,
              IWORK,
              IERR);
        }
      }
    }

    // Undo balancing on VL and VR and normalization
    // (Workspace: none needed)

    if (ILVL) {
      dggbak(BALANC, 'L', N, ILO.value, IHI.value, LSCALE, RSCALE, N, VL, LDVL,
          IERR);

      for (JC = 1; JC <= N; JC++) {
        if (ALPHAI[JC] < ZERO) continue;
        TEMP = ZERO;
        if (ALPHAI[JC] == ZERO) {
          for (JR = 1; JR <= N; JR++) {
            TEMP = max(TEMP, (VL[JR][JC].abs()));
          }
        } else {
          for (JR = 1; JR <= N; JR++) {
            TEMP = max(TEMP, VL[JR][JC].abs() + VL[JR][JC + 1].abs());
          }
        }
        if (TEMP < SMLNUM) continue;
        TEMP = ONE / TEMP;
        if (ALPHAI[JC] == ZERO) {
          for (JR = 1; JR <= N; JR++) {
            VL[JR][JC] *= TEMP;
          }
        } else {
          for (JR = 1; JR <= N; JR++) {
            VL[JR][JC] *= TEMP;
            VL[JR][JC + 1] *= TEMP;
          }
        }
      }
    }
    if (ILVR) {
      dggbak(BALANC, 'R', N, ILO.value, IHI.value, LSCALE, RSCALE, N, VR, LDVR,
          IERR);
      for (JC = 1; JC <= N; JC++) {
        if (ALPHAI[JC] < ZERO) continue;
        TEMP = ZERO;
        if (ALPHAI[JC] == ZERO) {
          for (JR = 1; JR <= N; JR++) {
            TEMP = max(TEMP, VR[JR][JC].abs());
          }
        } else {
          for (JR = 1; JR <= N; JR++) {
            TEMP = max(TEMP, VR[JR][JC].abs() + VR[JR][JC + 1].abs());
          }
        }
        if (TEMP < SMLNUM) continue;
        TEMP = ONE / TEMP;
        if (ALPHAI[JC] == ZERO) {
          for (JR = 1; JR <= N; JR++) {
            VR[JR][JC] *= TEMP;
          }
        } else {
          for (JR = 1; JR <= N; JR++) {
            VR[JR][JC] *= TEMP;
            VR[JR][JC + 1] *= TEMP;
          }
        }
      }
    }
  } finally {
    // Undo scaling if necessary

    if (ILASCL) {
      dlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAR.asMatrix(N), N, IERR);
      dlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAI.asMatrix(N), N, IERR);
    }

    if (ILBSCL) {
      dlascl('G', 0, 0, BNRMTO, BNRM, N, 1, BETA.asMatrix(N), N, IERR);
    }

    WORK[1] = MAXWRK.toDouble();
  }
}
