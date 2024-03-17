import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/variants/qr/ll/zgeqrf.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zggbak.dart';
import 'package:lapack/src/zggbal.dart';
import 'package:lapack/src/zgghrd.dart';
import 'package:lapack/src/zhgeqz.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlascl.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/ztgevc.dart';
import 'package:lapack/src/ztgsna.dart';
import 'package:lapack/src/zungqr.dart';
import 'package:lapack/src/zunmqr.dart';

void zggevx(
  final String BALANC,
  final String JOBVL,
  final String JOBVR,
  final String SENSE,
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
  final Box<int> ILO,
  final Box<int> IHI,
  final Array<double> LSCALE_,
  final Array<double> RSCALE_,
  final Box<double> ABNRM,
  final Box<double> BBNRM,
  final Array<double> RCONDE_,
  final Array<double> RCONDV_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Array<bool> BWORK_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final VL = VL_.having(ld: LDVL);
  final VR = VR_.having(ld: LDVR);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  final BWORK = BWORK_.having();
  final ALPHA = ALPHA_.having();
  final BETA = BETA_.having();
  final LSCALE = LSCALE_.having();
  final RSCALE = RSCALE_.having();
  final RCONDE = RCONDE_.having();
  final RCONDV = RCONDV_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool ILASCL,
      ILBSCL,
      ILV,
      ILVL,
      ILVR,
      LQUERY,
      NOSCL,
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
      MINWRK;
  double ANRM, ANRMTO = 0, BIGNUM, BNRM, BNRMTO = 0, EPS, SMLNUM, TEMP;
  final LDUMMA = Array<bool>(1);
  final IERR = Box(0), IN = Box(0), M = Box(0);

  double ABS1(Complex X) => X.toDouble().abs() + X.imaginary.abs();

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
  if (!(NOSCL || lsame(BALANC, 'S') || lsame(BALANC, 'B'))) {
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
    INFO.value = -13;
  } else if (LDVR < 1 || (ILVR && LDVR < N)) {
    INFO.value = -15;
  }

  // Compute workspace
  //  (Note: Comments in the code beginning "Workspace:" describe the
  //   minimal amount of workspace needed at that point in the code,
  //   as well as the preferred amount for good performance.
  //   NB refers to the optimal block size for the immediately
  //   following subroutine, as returned by ILAENV. The workspace is
  //   computed assuming ILO.value = 1 and IHI.value = N, the worst case.)

  if (INFO.value == 0) {
    if (N == 0) {
      MINWRK = 1;
      MAXWRK = 1;
    } else {
      MINWRK = 2 * N;
      if (WANTSE) {
        MINWRK = 4 * N;
      } else if (WANTSV || WANTSB) {
        MINWRK = 2 * N * (N + 1);
      }
      MAXWRK = MINWRK;
      MAXWRK = max(MAXWRK, N + N * ilaenv(1, 'ZGEQRF', ' ', N, 1, N, 0));
      MAXWRK = max(MAXWRK, N + N * ilaenv(1, 'ZUNMQR', ' ', N, 1, N, 0));
      if (ILVL) {
        MAXWRK = max(MAXWRK, N + N * ilaenv(1, 'ZUNGQR', ' ', N, 1, N, 0));
      }
    }
    WORK[1] = MAXWRK.toComplex();

    if (LWORK < MINWRK && !LQUERY) {
      INFO.value = -25;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZGGEVX', -INFO.value);
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

  // Permute and/or balance the matrix pair (A,B)
  // (Real Workspace: need 6*N if BALANC = 'S' or 'B', 1 otherwise)

  zggbal(BALANC, N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, RWORK, IERR);

  // Compute ABNRM.value and BBNRM.value

  ABNRM.value = zlange('1', N, N, A, LDA, RWORK(1));
  if (ILASCL) {
    RWORK[1] = ABNRM.value;
    dlascl('G', 0, 0, ANRMTO, ANRM, 1, 1, RWORK(1).asMatrix(1), 1, IERR);
    ABNRM.value = RWORK[1];
  }

  BBNRM.value = zlange('1', N, N, B, LDB, RWORK(1));
  if (ILBSCL) {
    RWORK[1] = BBNRM.value;
    dlascl('G', 0, 0, BNRMTO, BNRM, 1, 1, RWORK(1).asMatrix(1), 1, IERR);
    BBNRM.value = RWORK[1];
  }

  // Reduce B to triangular form (QR decomposition of B)
  // (Complex Workspace: need N, prefer N*NB )

  IROWS = IHI.value + 1 - ILO.value;
  if (ILV || !WANTSN) {
    ICOLS = N + 1 - ILO.value;
  } else {
    ICOLS = IROWS;
  }
  ITAU = 1;
  IWRK = ITAU + IROWS;
  zgeqrf(IROWS, ICOLS, B(ILO.value, ILO.value), LDB, WORK(ITAU), WORK(IWRK),
      LWORK + 1 - IWRK, IERR);

  // Apply the unitary transformation to A
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

  // Initialize VL and/or VR
  // (Workspace: need N, prefer N*NB)

  if (ILVL) {
    zlaset('Full', N, N, Complex.zero, Complex.one, VL, LDVL);
    if (IROWS > 1) {
      zlacpy('L', IROWS - 1, IROWS - 1, B(ILO.value + 1, ILO.value), LDB,
          VL(ILO.value + 1, ILO.value), LDVL);
    }
    zungqr(IROWS, IROWS, IROWS, VL(ILO.value, ILO.value), LDVL, WORK(ITAU),
        WORK(IWRK), LWORK + 1 - IWRK, IERR);
  }

  if (ILVR) zlaset('Full', N, N, Complex.zero, Complex.one, VR, LDVR);

  // Reduce to generalized Hessenberg form
  // (Workspace: none needed)

  if (ILV || !WANTSN) {
    // Eigenvectors requested -- work on whole matrix.

    zgghrd(JOBVL, JOBVR, N, ILO.value, IHI.value, A, LDA, B, LDB, VL, LDVL, VR,
        LDVR, IERR);
  } else {
    zgghrd('N', 'N', IROWS, 1, IROWS, A(ILO.value, ILO.value), LDA,
        B(ILO.value, ILO.value), LDB, VL, LDVL, VR, LDVR, IERR);
  }

  // Perform QZ algorithm (Compute eigenvalues, and optionally, the
  // Schur forms and Schur vectors)
  // (Complex Workspace: need N)
  // (Real Workspace: need N)

  IWRK = ITAU;
  if (ILV || !WANTSN) {
    CHTEMP = 'S';
  } else {
    CHTEMP = 'E';
  }

  failure:
  while (true) {
    zhgeqz(CHTEMP, JOBVL, JOBVR, N, ILO.value, IHI.value, A, LDA, B, LDB, ALPHA,
        BETA, VL, LDVL, VR, LDVR, WORK(IWRK), LWORK + 1 - IWRK, RWORK, IERR);
    if (IERR.value != 0) {
      if (IERR.value > 0 && IERR.value <= N) {
        INFO.value = IERR.value;
      } else if (IERR.value > N && IERR.value <= 2 * N) {
        INFO.value = IERR.value - N;
      } else {
        INFO.value = N + 1;
      }
      break failure;
    }

    // Compute Eigenvectors and estimate condition numbers if desired
    // ZTGEVC: (Complex Workspace: need 2*N )
    //         (Real Workspace:    need 2*N )
    // ZTGSNA: (Complex Workspace: need 2*N*N if SENSE='V' or 'B')
    //         (Integer Workspace: need N+2 )

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

        ztgevc(CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, N,
            IN, WORK(IWRK), RWORK, IERR);
        if (IERR.value != 0) {
          INFO.value = N + 2;
          break failure;
        }
      }

      if (!WANTSN) {
        // compute eigenvectors (ZTGEVC) and estimate condition
        // numbers (ZTGSNA). Note that the definition of the condition
        // number is not invariant under transformation (u,v) to
        // (Q*u, Z*v), where (u,v) are eigenvectors of the generalized
        // Schur form (S,T), Q and Z are orthogonal matrices. In order
        // to avoid using extra 2*N*N workspace, we have to
        // re-calculate eigenvectors and estimate the condition numbers
        // one at a time.

        for (I = 1; I <= N; I++) {
          for (J = 1; J <= N; J++) {
            BWORK[J] = false;
          }
          BWORK[I] = true;

          IWRK = N + 1;
          IWRK1 = IWRK + N;

          if (WANTSE || WANTSB) {
            ztgevc('B', 'S', BWORK, N, A, LDA, B, LDB, WORK(1).asMatrix(N), N,
                WORK(IWRK).asMatrix(N), N, 1, M, WORK(IWRK1), RWORK, IERR);
            if (IERR.value != 0) {
              INFO.value = N + 2;
              break failure;
            }
          }

          ztgsna(
              SENSE,
              'S',
              BWORK,
              N,
              A,
              LDA,
              B,
              LDB,
              WORK(1).asMatrix(N),
              N,
              WORK(IWRK).asMatrix(N),
              N,
              RCONDE(I),
              RCONDV(I),
              1,
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
      zggbak(BALANC, 'L', N, ILO.value, IHI.value, LSCALE, RSCALE, N, VL, LDVL,
          IERR);

      for (JC = 1; JC <= N; JC++) {
        TEMP = ZERO;
        for (JR = 1; JR <= N; JR++) {
          TEMP = max(TEMP, ABS1(VL[JR][JC]));
        }
        if (TEMP < SMLNUM) continue;
        TEMP = ONE / TEMP;
        for (JR = 1; JR <= N; JR++) {
          VL[JR][JC] *= TEMP.toComplex();
        }
      }
    }

    if (ILVR) {
      zggbak(BALANC, 'R', N, ILO.value, IHI.value, LSCALE, RSCALE, N, VR, LDVR,
          IERR);
      for (JC = 1; JC <= N; JC++) {
        TEMP = ZERO;
        for (JR = 1; JR <= N; JR++) {
          TEMP = max(TEMP, ABS1(VR[JR][JC]));
        }
        if (TEMP < SMLNUM) continue;
        TEMP = ONE / TEMP;
        for (JR = 1; JR <= N; JR++) {
          VR[JR][JC] *= TEMP.toComplex();
        }
      }
    }
  }

  // Undo scaling if necessary

  if (ILASCL) zlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA.asMatrix(N), N, IERR);

  if (ILBSCL) zlascl('G', 0, 0, BNRMTO, BNRM, N, 1, BETA.asMatrix(N), N, IERR);

  WORK[1] = MAXWRK.toComplex();
}
