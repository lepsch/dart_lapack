import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/lsame.dart';
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

void dgeesx(
  final String JOBVS,
  final String SORT,
  final bool Function(double, double) SELECT,
  final String SENSE,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Box<int> SDIM,
  final Array<double> WR_,
  final Array<double> WI_,
  final Matrix<double> VS_,
  final int LDVS,
  final Box<double> RCONDE,
  final Box<double> RCONDV,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Array<bool> BWORK_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final WR = WR_.dim();
  final WI = WI_.dim();
  final VS = VS_.dim(LDVS);
  final WORK = WORK_.dim();
  final IWORK = IWORK_.dim();
  final BWORK = BWORK_.dim();
  const ZERO = 0.0, ONE = 1.0;
  bool CURSL,
      LASTSL,
      LQUERY,
      LST2SL,
      SCALEA,
      WANTSB,
      WANTSE,
      WANTSN,
      WANTST,
      WANTSV,
      WANTVS;
  int HSWORK,
      I,
      I1,
      I2,
      IBAL,
      INXT,
      IP,
      ITAU,
      IWRK,
      LIWRK,
      LWRK,
      MAXWRK = 0,
      MINWRK;
  double ANRM, BIGNUM, CSCALE = 0, EPS, SMLNUM;
  final DUM = Array<double>(1);
  final IEVAL = Box(0),
      IERR = Box(0),
      IHI = Box(0),
      ILO = Box(0),
      ICOND = Box(0);

  // Test the input arguments

  INFO.value = 0;
  WANTVS = lsame(JOBVS, 'V');
  WANTST = lsame(SORT, 'S');
  WANTSN = lsame(SENSE, 'N');
  WANTSE = lsame(SENSE, 'E');
  WANTSV = lsame(SENSE, 'V');
  WANTSB = lsame(SENSE, 'B');
  LQUERY = (LWORK == -1 || LIWORK == -1);

  if ((!WANTVS) && (!lsame(JOBVS, 'N'))) {
    INFO.value = -1;
  } else if ((!WANTST) && (!lsame(SORT, 'N'))) {
    INFO.value = -2;
  } else if (!(WANTSN || WANTSE || WANTSV || WANTSB) || (!WANTST && !WANTSN)) {
    INFO.value = -4;
  } else if (N < 0) {
    INFO.value = -5;
  } else if (LDA < max(1, N)) {
    INFO.value = -7;
  } else if (LDVS < 1 || (WANTVS && LDVS < N)) {
    INFO.value = -12;
  }

  // Compute workspace
  //  (Note: Comments in the code beginning "RWorkspace:" describe the
  //   minimal amount of real workspace needed at that point in the
  //   code, as well as the preferred amount for good performance.
  //   IWorkspace refers to integer workspace.
  //   NB refers to the optimal block size for the immediately
  //   following subroutine, as returned by ILAENV.
  //   HSWORK refers to the workspace preferred by DHSEQR, as
  //   calculated below. HSWORK is computed assuming ILO.value=1 and IHI.value=N,
  //   the worst case.
  //   If SENSE = 'E', 'V' or 'B', then the amount of workspace needed
  //   depends on SDIM.value, which is computed by the routine DTRSEN later
  //   in the code.)

  if (INFO.value == 0) {
    LIWRK = 1;
    if (N == 0) {
      MINWRK = 1;
      LWRK = 1;
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
      LWRK = MAXWRK;
      if (!WANTSN) LWRK = max(LWRK, N + (N * N) ~/ 2);
      if (WANTSV || WANTSB) LIWRK = (N * N) ~/ 4;
    }
    IWORK[1] = LIWRK;
    WORK[1] = LWRK.toDouble();

    if (LWORK < MINWRK && !LQUERY) {
      INFO.value = -16;
    } else if (LIWORK < 1 && !LQUERY) {
      INFO.value = -18;
    }
  }

  if (INFO.value != 0) {
    xerbla('DGEESX', -INFO.value);
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
  // (RWorkspace: need N)

  IBAL = 1;
  dgebal('P', N, A, LDA, ILO, IHI, WORK(IBAL), IERR);

  // Reduce to upper Hessenberg form
  // (RWorkspace: need 3*N, prefer 2*N+N*NB)

  ITAU = N + IBAL;
  IWRK = N + ITAU;
  dgehrd(N, ILO.value, IHI.value, A, LDA, WORK(ITAU), WORK(IWRK),
      LWORK - IWRK + 1, IERR);

  if (WANTVS) {
    // Copy Householder vectors to VS

    dlacpy('L', N, N, A, LDA, VS, LDVS);

    // Generate orthogonal matrix in VS
    // (RWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)

    dorghr(N, ILO.value, IHI.value, VS, LDVS, WORK(ITAU), WORK(IWRK),
        LWORK - IWRK + 1, IERR);
  }

  SDIM.value = 0;

  // Perform QR iteration, accumulating Schur vectors in VS if desired
  // (RWorkspace: need N+1, prefer N+HSWORK (see comments) )

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

    // Reorder eigenvalues, transform Schur vectors, and compute
    // reciprocal condition numbers
    // (RWorkspace: if SENSE is not 'N', need N+2*SDIM.value*(N-SDIM.value)
    //              otherwise, need N )
    // (IWorkspace: if SENSE is 'V' or 'B', need SDIM.value*(N-SDIM.value)
    //              otherwise, need 0 )

    dtrsen(SENSE, JOBVS, BWORK, N, A, LDA, VS, LDVS, WR, WI, SDIM, RCONDE,
        RCONDV, WORK(IWRK), LWORK - IWRK + 1, IWORK, LIWORK, ICOND);
    if (!WANTSN) MAXWRK = max(MAXWRK, N + 2 * SDIM.value * (N - SDIM.value));
    if (ICOND.value == -15) {
      // Not enough real workspace

      INFO.value = -16;
    } else if (ICOND.value == -17) {
      // Not enough integer workspace

      INFO.value = -18;
    } else if (ICOND.value > 0) {
      // DTRSEN failed to reorder or to restore standard Schur form

      INFO.value = ICOND.value + N;
    }
  }

  if (WANTVS) {
    // Undo balancing
    // (RWorkspace: need N)

    dgebak('P', 'R', N, ILO.value, IHI.value, WORK(IBAL), N, VS, LDVS, IERR);
  }

  if (SCALEA) {
    // Undo scaling for the Schur form of A

    dlascl('H', 0, 0, CSCALE, ANRM, N, N, A, LDA, IERR);
    dcopy(N, A.asArray(), LDA + 1, WR, 1);
    if ((WANTSV || WANTSB) && INFO.value == 0) {
      DUM[1] = RCONDV.value;
      dlascl('G', 0, 0, CSCALE, ANRM, 1, 1, DUM.asMatrix(1), 1, IERR);
      RCONDV.value = DUM[1];
    }
    if (CSCALE == SMLNUM) {
      // If scaling back towards underflow, adjust WI if an
      // offdiagonal element of a 2-by-2 block in the Schur form
      // underflows.

      if (IEVAL.value > 0) {
        I1 = IEVAL.value + 1;
        I2 = IHI.value - 1;
        dlascl(
            'G', 0, 0, CSCALE, ANRM, ILO.value - 1, 1, WI.asMatrix(N), N, IERR);
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
            if (I > 1) {
              dswap(I - 1, A(1, I).asArray(), 1, A(1, I + 1).asArray(), 1);
            }
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
    final LDWI = max(N - IEVAL.value, 1);
    dlascl('G', 0, 0, CSCALE, ANRM, N - IEVAL.value, 1,
        WI(IEVAL.value + 1).asMatrix(LDWI), LDWI, IERR);
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
        if (CURSL) SDIM.value = SDIM.value + 1;
        IP = 0;
        if (CURSL && !LASTSL) INFO.value = N + 2;
      } else {
        if (IP == 1) {
          // Last eigenvalue of conjugate pair

          CURSL = CURSL || LASTSL;
          LASTSL = CURSL;
          if (CURSL) SDIM.value = SDIM.value + 2;
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
  if (WANTSV || WANTSB) {
    IWORK[1] = max(1, SDIM.value * (N - SDIM.value));
  } else {
    IWORK[1] = 1;
  }
}
