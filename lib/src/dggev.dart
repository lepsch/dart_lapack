import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgeqrf.dart';
import 'package:lapack/src/dggbak.dart';
import 'package:lapack/src/dggbal.dart';
import 'package:lapack/src/dgghrd.dart';
import 'package:lapack/src/dhgeqz.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dorgqr.dart';
import 'package:lapack/src/dormqr.dart';
import 'package:lapack/src/dtgevc.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dggev(
  final String JOBVL,
  final String JOBVR,
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
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final B = B_.dim(LDB);
  final ALPHAR = ALPHAR_.dim();
  final ALPHAI = ALPHAI_.dim();
  final BETA = BETA_.dim();
  final VL = VL_.dim(LDVL);
  final VR = VR_.dim(LDVR);
  final WORK = WORK_.dim();
  const ZERO = 0.0, ONE = 1.0;
  bool ILASCL, ILBSCL, ILV, ILVL, ILVR, LQUERY;
  String CHTEMP;
  int ICOLS,
      IJOBVL,
      IJOBVR,
      ILEFT,
      IRIGHT,
      IROWS,
      ITAU,
      IWRK,
      JC,
      JR,
      MAXWRK = 0,
      MINWRK;
  double ANRM, ANRMTO = 0, BIGNUM, BNRM, BNRMTO = 0, EPS, SMLNUM, TEMP;
  final LDUMMA = Array<bool>(1);
  final IERR = Box(0), ILO = Box(0), IHI = Box(0), IN = Box(0);

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
    INFO.value = -12;
  } else if (LDVR < 1 || (ILVR && LDVR < N)) {
    INFO.value = -14;
  }

  // Compute workspace
  // (Note: Comments in the code beginning "Workspace:" describe the
  // minimal amount of workspace needed at that point in the code,
  // as well as the preferred amount for good performance.
  // NB refers to the optimal block size for the immediately
  // following subroutine, as returned by ILAENV. The workspace is
  // computed assuming ILO.value = 1 and IHI.value = N, the worst case.)

  if (INFO.value == 0) {
    MINWRK = max(1, 8 * N);
    MAXWRK = max(1, N * (7 + ilaenv(1, 'DGEQRF', ' ', N, 1, N, 0)));
    MAXWRK = max(MAXWRK, N * (7 + ilaenv(1, 'DORMQR', ' ', N, 1, N, 0)));
    if (ILVL) {
      MAXWRK = max(MAXWRK, N * (7 + ilaenv(1, 'DORGQR', ' ', N, 1, N, -1)));
    }
    WORK[1] = MAXWRK.toDouble();

    if (LWORK < MINWRK && !LQUERY) INFO.value = -16;
  }

  if (INFO.value != 0) {
    xerbla('DGGEV ', -INFO.value);
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

  // Permute the matrices A, B to isolate eigenvalues if possible
  // (Workspace: need 6*N)

  ILEFT = 1;
  IRIGHT = N + 1;
  IWRK = IRIGHT + N;
  dggbal('P', N, A, LDA, B, LDB, ILO, IHI, WORK(ILEFT), WORK(IRIGHT),
      WORK(IWRK), IERR);

  // Reduce B to triangular form (QR decomposition of B)
  // (Workspace: need N, prefer N*NB)

  IROWS = IHI.value + 1 - ILO.value;
  if (ILV) {
    ICOLS = N + 1 - ILO.value;
  } else {
    ICOLS = IROWS;
  }
  ITAU = IWRK;
  IWRK = ITAU + IROWS;
  dgeqrf(IROWS, ICOLS, B(ILO.value, ILO.value), LDB, WORK(ITAU), WORK(IWRK),
      LWORK + 1 - IWRK, IERR);

  // Apply the orthogonal transformation to matrix A
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

  // Initialize VL
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

  // Initialize VR

  if (ILVR) dlaset('Full', N, N, ZERO, ONE, VR, LDVR);

  // Reduce to generalized Hessenberg form
  // (Workspace: none needed)

  if (ILV) {
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

  IWRK = ITAU;
  if (ILV) {
    CHTEMP = 'S';
  } else {
    CHTEMP = 'E';
  }
  dhgeqz(CHTEMP, JOBVL, JOBVR, N, ILO.value, IHI.value, A, LDA, B, LDB, ALPHAR,
      ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK(IWRK), LWORK + 1 - IWRK, IERR);
  if (IERR.value != 0) {
    if (IERR.value > 0 && IERR.value <= N) {
      INFO.value = IERR.value;
    } else if (IERR.value > N && IERR.value <= 2 * N) {
      INFO.value = IERR.value - N;
    } else {
      INFO.value = N + 1;
    }
  } else {
    // Compute Eigenvectors
    // (Workspace: need 6*N)

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
      dtgevc(CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, N, IN,
          WORK(IWRK), IERR);
      if (IERR.value != 0) {
        INFO.value = N + 2;
      } else {
        // Undo balancing on VL and VR and normalization
        // (Workspace: none needed)

        if (ILVL) {
          dggbak('P', 'L', N, ILO.value, IHI.value, WORK(ILEFT), WORK(IRIGHT),
              N, VL, LDVL, IERR);
          for (JC = 1; JC <= N; JC++) {
            if (ALPHAI[JC] < ZERO) continue;
            TEMP = ZERO;
            if (ALPHAI[JC] == ZERO) {
              for (JR = 1; JR <= N; JR++) {
                TEMP = max(TEMP, (VL[JR][JC]).abs());
              }
            } else {
              for (JR = 1; JR <= N; JR++) {
                TEMP = max(TEMP, (VL[JR][JC]).abs() + (VL[JR][JC + 1]).abs());
              }
            }
            if (TEMP < SMLNUM) continue;
            TEMP = ONE / TEMP;
            if (ALPHAI[JC] == ZERO) {
              for (JR = 1; JR <= N; JR++) {
                VL[JR][JC] = VL[JR][JC] * TEMP;
              }
            } else {
              for (JR = 1; JR <= N; JR++) {
                VL[JR][JC] = VL[JR][JC] * TEMP;
                VL[JR][JC + 1] = VL[JR][JC + 1] * TEMP;
              }
            }
          }
        }
        if (ILVR) {
          dggbak('P', 'R', N, ILO.value, IHI.value, WORK(ILEFT), WORK(IRIGHT),
              N, VR, LDVR, IERR);
          for (JC = 1; JC <= N; JC++) {
            if (ALPHAI[JC] < ZERO) continue;
            TEMP = ZERO;
            if (ALPHAI[JC] == ZERO) {
              for (JR = 1; JR <= N; JR++) {
                TEMP = max(TEMP, (VR[JR][JC]).abs());
              }
            } else {
              for (JR = 1; JR <= N; JR++) {
                TEMP = max(TEMP, (VR[JR][JC]).abs() + (VR[JR][JC + 1]).abs());
              }
            }
            if (TEMP < SMLNUM) continue;
            TEMP = ONE / TEMP;
            if (ALPHAI[JC] == ZERO) {
              for (JR = 1; JR <= N; JR++) {
                VR[JR][JC] = VR[JR][JC] * TEMP;
              }
            } else {
              for (JR = 1; JR <= N; JR++) {
                VR[JR][JC] = VR[JR][JC] * TEMP;
                VR[JR][JC + 1] = VR[JR][JC + 1] * TEMP;
              }
            }
          }
        }
      }

      // End of eigenvector calculation
    }
  }

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
