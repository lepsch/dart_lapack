import 'dart:math';

import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/blas/drot.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/blas/lsame.dart';
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
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgeev(
  final String JOBVL,
  final String JOBVR,
  final int N,
  final Matrix<double> A,
  final int LDA,
  final Array<double> WR,
  final Array<double> WI,
  final Matrix<double> VL,
  final int LDVL,
  final Matrix<double> VR,
  final int LDVR,
  final Array<double> WORK,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  bool LQUERY, SCALEA, WANTVL, WANTVR;
  String SIDE = '';
  int HSWORK, I, IBAL, ITAU, IWRK, K, LWORK_TREVC, MAXWRK = 0, MINWRK;
  double ANRM, BIGNUM, CSCALE = 0, EPS, SCL, SMLNUM;
  final SELECT = Array<bool>(1);
  final DUM = Array<double>(1);
  final NOUT = Box(0), IERR = Box(0), IHI = Box(0), ILO = Box(0);
  final CS = Box(0.0), R = Box(0.0), SN = Box(0.0);

  // Test the input arguments

  INFO.value = 0;
  LQUERY = (LWORK == -1);
  WANTVL = lsame(JOBVL, 'V');
  WANTVR = lsame(JOBVR, 'V');
  if ((!WANTVL) && (!lsame(JOBVL, 'N'))) {
    INFO.value = -1;
  } else if ((!WANTVR) && (!lsame(JOBVR, 'N'))) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDVL < 1 || (WANTVL && LDVL < N)) {
    INFO.value = -9;
  } else if (LDVR < 1 || (WANTVR && LDVR < N)) {
    INFO.value = -11;
  }

  // Compute workspace
  // (Note: Comments in the code beginning "Workspace:" describe the
  // minimal amount of workspace needed at that point in the code,
  // as well as the preferred amount for good performance.
  // NB refers to the optimal block size for the immediately
  // following subroutine, as returned by ILAENV.
  // HSWORK refers to the workspace preferred by DHSEQR, as
  // calculated below. HSWORK is computed assuming ILO.value=1 and IHI.value=N,
  // the worst case.)

  if (INFO.value == 0) {
    if (N == 0) {
      MINWRK = 1;
      MAXWRK = 1;
    } else {
      MAXWRK = 2 * N + N * ilaenv(1, 'DGEHRD', ' ', N, 1, N, 0);
      if (WANTVL) {
        MINWRK = 4 * N;
        MAXWRK = max(
          MAXWRK,
          2 * N + (N - 1) * ilaenv(1, 'DORGHR', ' ', N, 1, N, -1),
        );
        dhseqr('S', 'V', N, 1, N, A, LDA, WR, WI, VL, LDVL, WORK, -1, INFO);
        HSWORK = WORK[1].toInt();
        MAXWRK = max(MAXWRK, max(N + 1, N + HSWORK));
        dtrevc3(
          'L',
          'B',
          SELECT,
          N,
          A,
          LDA,
          VL,
          LDVL,
          VR,
          LDVR,
          N,
          NOUT,
          WORK,
          -1,
          IERR,
        );
        LWORK_TREVC = WORK[1].toInt();
        MAXWRK = max(MAXWRK, N + LWORK_TREVC);
        MAXWRK = max(MAXWRK, 4 * N);
      } else if (WANTVR) {
        MINWRK = 4 * N;
        MAXWRK = max(
          MAXWRK,
          2 * N + (N - 1) * ilaenv(1, 'DORGHR', ' ', N, 1, N, -1),
        );
        dhseqr('S', 'V', N, 1, N, A, LDA, WR, WI, VR, LDVR, WORK, -1, INFO);
        HSWORK = WORK[1].toInt();
        MAXWRK = max(MAXWRK, max(N + 1, N + HSWORK));
        dtrevc3(
          'R',
          'B',
          SELECT,
          N,
          A,
          LDA,
          VL,
          LDVL,
          VR,
          LDVR,
          N,
          NOUT,
          WORK,
          -1,
          IERR,
        );
        LWORK_TREVC = WORK[1].toInt();
        MAXWRK = max(MAXWRK, N + LWORK_TREVC);
        MAXWRK = max(MAXWRK, 4 * N);
      } else {
        MINWRK = 3 * N;
        dhseqr('E', 'N', N, 1, N, A, LDA, WR, WI, VR, LDVR, WORK, -1, INFO);
        HSWORK = WORK[1].toInt();
        MAXWRK = max(MAXWRK, max(N + 1, N + HSWORK));
      }
      MAXWRK = max(MAXWRK, MINWRK);
    }
    WORK[1] = MAXWRK.toDouble();

    if (LWORK < MINWRK && !LQUERY) {
      INFO.value = -13;
    }
  }

  if (INFO.value != 0) {
    xerbla('DGEEV ', -INFO.value);
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

  // Balance the matrix
  // (Workspace: need N)

  IBAL = 1;
  dgebal('B', N, A, LDA, ILO, IHI, WORK(IBAL), IERR);

  // Reduce to upper Hessenberg form
  // (Workspace: need 3*N, prefer 2*N+N*NB)

  ITAU = IBAL + N;
  IWRK = ITAU + N;
  dgehrd(
    N,
    ILO.value,
    IHI.value,
    A,
    LDA,
    WORK(ITAU),
    WORK(IWRK),
    LWORK - IWRK + 1,
    IERR,
  );

  if (WANTVL) {
    // Want left eigenvectors
    // Copy Householder vectors to VL

    SIDE = 'L';
    dlacpy('L', N, N, A, LDA, VL, LDVL);

    // Generate orthogonal matrix in VL
    // (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)

    dorghr(
      N,
      ILO.value,
      IHI.value,
      VL,
      LDVL,
      WORK(ITAU),
      WORK(IWRK),
      LWORK - IWRK + 1,
      IERR,
    );

    // Perform QR iteration, accumulating Schur vectors in VL
    // (Workspace: need N+1, prefer N+HSWORK (see comments) )

    IWRK = ITAU;
    dhseqr(
      'S',
      'V',
      N,
      ILO.value,
      IHI.value,
      A,
      LDA,
      WR,
      WI,
      VL,
      LDVL,
      WORK(IWRK),
      LWORK - IWRK + 1,
      INFO,
    );

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
    // (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)

    dorghr(
      N,
      ILO.value,
      IHI.value,
      VR,
      LDVR,
      WORK(ITAU),
      WORK(IWRK),
      LWORK - IWRK + 1,
      IERR,
    );

    // Perform QR iteration, accumulating Schur vectors in VR
    // (Workspace: need N+1, prefer N+HSWORK (see comments) )

    IWRK = ITAU;
    dhseqr(
      'S',
      'V',
      N,
      ILO.value,
      IHI.value,
      A,
      LDA,
      WR,
      WI,
      VR,
      LDVR,
      WORK(IWRK),
      LWORK - IWRK + 1,
      INFO,
    );
  } else {
    // Compute eigenvalues only
    // (Workspace: need N+1, prefer N+HSWORK (see comments) )

    IWRK = ITAU;
    dhseqr(
      'E',
      'N',
      N,
      ILO.value,
      IHI.value,
      A,
      LDA,
      WR,
      WI,
      VR,
      LDVR,
      WORK(IWRK),
      LWORK - IWRK + 1,
      INFO,
    );
  }

  // If INFO.value != 0 from DHSEQR, then quit

  if (INFO.value == 0) {
    if (WANTVL || WANTVR) {
      // Compute left and/or right eigenvectors
      // (Workspace: need 4*N, prefer N + N + 2*N*NB)

      dtrevc3(
        SIDE,
        'B',
        SELECT,
        N,
        A,
        LDA,
        VL,
        LDVL,
        VR,
        LDVR,
        N,
        NOUT,
        WORK(IWRK),
        LWORK - IWRK + 1,
        IERR,
      );
    }

    if (WANTVL) {
      // Undo balancing of left eigenvectors
      // (Workspace: need N)

      dgebak('B', 'L', N, ILO.value, IHI.value, WORK(IBAL), N, VL, LDVL, IERR);

      // Normalize left eigenvectors and make largest component real

      for (I = 1; I <= N; I++) {
        // 20
        if (WI[I] == ZERO) {
          SCL = ONE / dnrm2(N, VL(1, I).asArray(), 1);
          dscal(N, SCL, VL(1, I).asArray(), 1);
        } else if (WI[I] > ZERO) {
          SCL = ONE /
              dlapy2(
                dnrm2(N, VL(1, I).asArray(), 1),
                dnrm2(N, VL(1, I + 1).asArray(), 1),
              );
          dscal(N, SCL, VL(1, I).asArray(), 1);
          dscal(N, SCL, VL(1, I + 1).asArray(), 1);
          for (K = 1; K <= N; K++) {
            // 10
            WORK[IWRK + K - 1] =
                pow(VL[K][I], 2).toDouble() + pow(VL[K][I + 1], 2);
          } // 10
          K = idamax(N, WORK(IWRK), 1);
          dlartg(VL[K][I], VL[K][I + 1], CS, SN, R);
          drot(
            N,
            VL(1, I).asArray(),
            1,
            VL(1, I + 1).asArray(),
            1,
            CS.value,
            SN.value,
          );
          VL[K][I + 1] = ZERO;
        }
      } // 20
    }

    if (WANTVR) {
      // Undo balancing of right eigenvectors
      // (Workspace: need N)

      dgebak(
        'B',
        'R',
        N,
        ILO.value,
        IHI.value,
        WORK(IBAL),
        N,
        VR,
        LDVR,
        IERR,
      );

      // Normalize right eigenvectors and make largest component real

      for (I = 1; I <= N; I++) {
        // 40
        if (WI[I] == ZERO) {
          SCL = ONE / dnrm2(N, VR(1, I).asArray(), 1);
          dscal(N, SCL, VR(1, I).asArray(), 1);
        } else if (WI[I] > ZERO) {
          SCL = ONE /
              dlapy2(
                dnrm2(N, VR(1, I).asArray(), 1),
                dnrm2(N, VR(1, I + 1).asArray(), 1),
              );
          dscal(N, SCL, VR(1, I).asArray(), 1);
          dscal(N, SCL, VR(1, I + 1).asArray(), 1);
          for (K = 1; K <= N; K++) {
            // 30
            WORK[IWRK + K - 1] =
                pow(VR[K][I], 2).toDouble() + pow(VR[K][I + 1], 2);
          } // 30
          K = idamax(N, WORK(IWRK), 1);
          dlartg(VR[K][I], VR[K][I + 1], CS, SN, R);
          drot(
            N,
            VR(1, I).asArray(),
            1,
            VR(1, I + 1).asArray(),
            1,
            CS.value,
            SN.value,
          );
          VR[K][I + 1] = ZERO;
        }
      } // 40
    }
  } // 50

  // Undo scaling if necessary

  if (SCALEA) {
    var ld = max(N - INFO.value, 1);
    dlascl(
      'G',
      0,
      0,
      CSCALE,
      ANRM,
      N - INFO.value,
      1,
      WR(INFO.value + 1).asMatrix(ld),
      ld,
      IERR,
    );
    dlascl(
      'G',
      0,
      0,
      CSCALE,
      ANRM,
      N - INFO.value,
      1,
      WI(INFO.value + 1).asMatrix(ld),
      ld,
      IERR,
    );
    if (INFO.value > 0) {
      dlascl(
        'G',
        0,
        0,
        CSCALE,
        ANRM,
        ILO.value - 1,
        1,
        WR.asMatrix(N),
        N,
        IERR,
      );
      dlascl(
        'G',
        0,
        0,
        CSCALE,
        ANRM,
        ILO.value - 1,
        1,
        WI.asMatrix(N),
        N,
        IERR,
      );
    }
  }

  WORK[1] = MAXWRK.toDouble();
}
