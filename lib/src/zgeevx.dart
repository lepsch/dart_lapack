import 'dart:math';

import 'package:lapack/src/blas/dznrm2.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgebak.dart';
import 'package:lapack/src/zgebal.dart';
import 'package:lapack/src/zgehrd.dart';
import 'package:lapack/src/zhseqr.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlascl.dart';
import 'package:lapack/src/ztrevc3.dart';
import 'package:lapack/src/ztrsna.dart';
import 'package:lapack/src/zunghr.dart';

void zgeevx(
  final String BALANC,
  final String JOBVL,
  final String JOBVR,
  final String SENSE,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> W_,
  final Matrix<Complex> VL_,
  final int LDVL,
  final Matrix<Complex> VR_,
  final int LDVR,
  final Box<int> ILO,
  final Box<int> IHI,
  final Array<double> SCALE_,
  final Box<double> ABNRM,
  final Array<double> RCONDE_,
  final Array<double> RCONDV_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final VL = VL_.having(ld: LDVL);
  final VR = VR_.having(ld: LDVR);
  final W = W_.having();
  final SCALE = SCALE_.having();
  final RCONDE = RCONDE_.having();
  final RCONDV = RCONDV_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool LQUERY, SCALEA, WANTVL, WANTVR, WNTSNB, WNTSNE, WNTSNN, WNTSNV;
  String JOB, SIDE = '';
  int HSWORK, I, ITAU, IWRK, K, LWORK_TREVC, MAXWRK = 0, MINWRK;
  double ANRM, BIGNUM, CSCALE = 0, EPS, SCL, SMLNUM;
  Complex TMP;
  final SELECT = Array<bool>(1);
  final DUM = Array<double>(1);
  final IERR = Box(0), NOUT = Box(0), ICOND = Box(0);

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
  } else if ((!WANTVL) && (!lsame(JOBVL, 'N'))) {
    INFO.value = -2;
  } else if ((!WANTVR) && (!lsame(JOBVR, 'N'))) {
    INFO.value = -3;
  } else if (!(WNTSNN || WNTSNE || WNTSNB || WNTSNV) ||
      ((WNTSNE || WNTSNB) && !(WANTVL && WANTVR))) {
    INFO.value = -4;
  } else if (N < 0) {
    INFO.value = -5;
  } else if (LDA < max(1, N)) {
    INFO.value = -7;
  } else if (LDVL < 1 || (WANTVL && LDVL < N)) {
    INFO.value = -10;
  } else if (LDVR < 1 || (WANTVR && LDVR < N)) {
    INFO.value = -12;
  }

  // Compute workspace
  //  (Note: Comments in the code beginning "Workspace:" describe the
  //   minimal amount of workspace needed at that point in the code,
  //   as well as the preferred amount for good performance.
  //   CWorkspace refers to complex workspace, and RWorkspace to real
  //   workspace. NB refers to the optimal block size for the
  //   immediately following subroutine, as returned by ILAENV.
  //   HSWORK refers to the workspace preferred by ZHSEQR, as
  //   calculated below. HSWORK is computed assuming ILO.value=1 and IHI.value=N,
  //   the worst case.)

  if (INFO.value == 0) {
    if (N == 0) {
      MINWRK = 1;
      MAXWRK = 1;
    } else {
      MAXWRK = N + N * ilaenv(1, 'ZGEHRD', ' ', N, 1, N, 0);

      if (WANTVL) {
        ztrevc3('L', 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK,
            -1, RWORK, -1, IERR);
        LWORK_TREVC = WORK[1].toInt();
        MAXWRK = max(MAXWRK, LWORK_TREVC);
        zhseqr('S', 'V', N, 1, N, A, LDA, W, VL, LDVL, WORK, -1, INFO);
      } else if (WANTVR) {
        ztrevc3('R', 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK,
            -1, RWORK, -1, IERR);
        LWORK_TREVC = WORK[1].toInt();
        MAXWRK = max(MAXWRK, LWORK_TREVC);
        zhseqr('S', 'V', N, 1, N, A, LDA, W, VR, LDVR, WORK, -1, INFO);
      } else {
        if (WNTSNN) {
          zhseqr('E', 'N', N, 1, N, A, LDA, W, VR, LDVR, WORK, -1, INFO);
        } else {
          zhseqr('S', 'N', N, 1, N, A, LDA, W, VR, LDVR, WORK, -1, INFO);
        }
      }
      HSWORK = WORK[1].toInt();

      if ((!WANTVL) && (!WANTVR)) {
        MINWRK = 2 * N;
        if (!(WNTSNN || WNTSNE)) MINWRK = max(MINWRK, N * N + 2 * N);
        MAXWRK = max(MAXWRK, HSWORK);
        if (!(WNTSNN || WNTSNE)) MAXWRK = max(MAXWRK, N * N + 2 * N);
      } else {
        MINWRK = 2 * N;
        if (!(WNTSNN || WNTSNE)) MINWRK = max(MINWRK, N * N + 2 * N);
        MAXWRK = max(MAXWRK, HSWORK);
        MAXWRK =
            max(MAXWRK, N + (N - 1) * ilaenv(1, 'ZUNGHR', ' ', N, 1, N, -1));
        if (!(WNTSNN || WNTSNE)) MAXWRK = max(MAXWRK, N * N + 2 * N);
        MAXWRK = max(MAXWRK, 2 * N);
      }
      MAXWRK = max(MAXWRK, MINWRK);
    }
    WORK[1] = MAXWRK.toComplex();

    if (LWORK < MINWRK && !LQUERY) {
      INFO.value = -20;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZGEEVX', -INFO.value);
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
  ANRM = zlange('M', N, N, A, LDA, DUM);
  SCALEA = false;
  if (ANRM > ZERO && ANRM < SMLNUM) {
    SCALEA = true;
    CSCALE = SMLNUM;
  } else if (ANRM > BIGNUM) {
    SCALEA = true;
    CSCALE = BIGNUM;
  }
  if (SCALEA) zlascl('G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR);

  // Balance the matrix and compute ABNRM.value

  zgebal(BALANC, N, A, LDA, ILO, IHI, SCALE, IERR);
  ABNRM.value = zlange('1', N, N, A, LDA, DUM);
  if (SCALEA) {
    DUM[1] = ABNRM.value;
    dlascl('G', 0, 0, CSCALE, ANRM, 1, 1, DUM.asMatrix(1), 1, IERR);
    ABNRM.value = DUM[1];
  }

  // Reduce to upper Hessenberg form
  // (CWorkspace: need 2*N, prefer N+N*NB)
  // (RWorkspace: none)

  ITAU = 1;
  IWRK = ITAU + N;
  zgehrd(N, ILO.value, IHI.value, A, LDA, WORK(ITAU), WORK(IWRK),
      LWORK - IWRK + 1, IERR);

  if (WANTVL) {
    // Want left eigenvectors
    // Copy Householder vectors to VL

    SIDE = 'L';
    zlacpy('L', N, N, A, LDA, VL, LDVL);

    // Generate unitary matrix in VL
    // (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
    // (RWorkspace: none)

    zunghr(N, ILO.value, IHI.value, VL, LDVL, WORK(ITAU), WORK(IWRK),
        LWORK - IWRK + 1, IERR);

    // Perform QR iteration, accumulating Schur vectors in VL
    // (CWorkspace: need 1, prefer HSWORK (see comments) )
    // (RWorkspace: none)

    IWRK = ITAU;
    zhseqr('S', 'V', N, ILO.value, IHI.value, A, LDA, W, VL, LDVL, WORK(IWRK),
        LWORK - IWRK + 1, INFO);

    if (WANTVR) {
      // Want left and right eigenvectors
      // Copy Schur vectors to VR

      SIDE = 'B';
      zlacpy('F', N, N, VL, LDVL, VR, LDVR);
    }
  } else if (WANTVR) {
    // Want right eigenvectors
    // Copy Householder vectors to VR

    SIDE = 'R';
    zlacpy('L', N, N, A, LDA, VR, LDVR);

    // Generate unitary matrix in VR
    // (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
    // (RWorkspace: none)

    zunghr(N, ILO.value, IHI.value, VR, LDVR, WORK(ITAU), WORK(IWRK),
        LWORK - IWRK + 1, IERR);

    // Perform QR iteration, accumulating Schur vectors in VR
    // (CWorkspace: need 1, prefer HSWORK (see comments) )
    // (RWorkspace: none)

    IWRK = ITAU;
    zhseqr('S', 'V', N, ILO.value, IHI.value, A, LDA, W, VR, LDVR, WORK(IWRK),
        LWORK - IWRK + 1, INFO);
  } else {
    // Compute eigenvalues only
    // If condition numbers desired, compute Schur form

    if (WNTSNN) {
      JOB = 'E';
    } else {
      JOB = 'S';
    }

    // (CWorkspace: need 1, prefer HSWORK (see comments) )
    // (RWorkspace: none)

    IWRK = ITAU;
    zhseqr(JOB, 'N', N, ILO.value, IHI.value, A, LDA, W, VR, LDVR, WORK(IWRK),
        LWORK - IWRK + 1, INFO);
  }

  // If INFO.value != 0 from ZHSEQR, then quit

  if (INFO.value == 0) {
    if (WANTVL || WANTVR) {
      // Compute left and/or right eigenvectors
      // (CWorkspace: need 2*N, prefer N + 2*N*NB)
      // (RWorkspace: need N)

      ztrevc3(SIDE, 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT,
          WORK(IWRK), LWORK - IWRK + 1, RWORK, N, IERR);
    }

    // Compute condition numbers if desired
    // (CWorkspace: need N*N+2*N unless SENSE = 'E')
    // (RWorkspace: need 2*N unless SENSE = 'E')

    if (!WNTSNN) {
      ztrsna(SENSE, 'A', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, RCONDE, RCONDV,
          N, NOUT, WORK(IWRK).asMatrix(N), N, RWORK, ICOND);
    }

    if (WANTVL) {
      // Undo balancing of left eigenvectors

      zgebak(BALANC, 'L', N, ILO.value, IHI.value, SCALE, N, VL, LDVL, IERR);

      // Normalize left eigenvectors and make largest component real

      for (I = 1; I <= N; I++) {
        SCL = ONE / dznrm2(N, VL(1, I).asArray(), 1);
        zdscal(N, SCL, VL(1, I).asArray(), 1);
        for (K = 1; K <= N; K++) {
          RWORK[K] = pow(VL[K][I].toDouble(), 2).toDouble() +
              pow(VL[K][I].imaginary, 2);
        }
        K = idamax(N, RWORK, 1);
        TMP = VL[K][I].conjugate() / sqrt(RWORK[K]).toComplex();
        zscal(N, TMP, VL(1, I).asArray(), 1);
        VL[K][I] = VL[K][I].toDouble().toComplex();
      }
    }

    if (WANTVR) {
      // Undo balancing of right eigenvectors

      zgebak(BALANC, 'R', N, ILO.value, IHI.value, SCALE, N, VR, LDVR, IERR);

      // Normalize right eigenvectors and make largest component real

      for (I = 1; I <= N; I++) {
        SCL = ONE / dznrm2(N, VR(1, I).asArray(), 1);
        zdscal(N, SCL, VR(1, I).asArray(), 1);
        for (K = 1; K <= N; K++) {
          RWORK[K] = pow(VR[K][I].toDouble(), 2).toDouble() +
              pow(VR[K][I].imaginary, 2);
        }
        K = idamax(N, RWORK, 1);
        TMP = VR[K][I].conjugate() / sqrt(RWORK[K]).toComplex();
        zscal(N, TMP, VR(1, I).asArray(), 1);
        VR[K][I] = VR[K][I].toDouble().toComplex();
      }
    }

    // Undo scaling if necessary
  }
  if (SCALEA) {
    zlascl(
        'G',
        0,
        0,
        CSCALE,
        ANRM,
        N - INFO.value,
        1,
        W(INFO.value + 1).asMatrix(max(N - INFO.value, 1)),
        max(N - INFO.value, 1),
        IERR);
    if (INFO.value == 0) {
      if ((WNTSNV || WNTSNB) && ICOND.value == 0) {
        dlascl('G', 0, 0, CSCALE, ANRM, N, 1, RCONDV.asMatrix(N), N, IERR);
      }
    } else {
      zlascl('G', 0, 0, CSCALE, ANRM, ILO.value - 1, 1, W.asMatrix(N), N, IERR);
    }
  }

  WORK[1] = MAXWRK.toComplex();
}
