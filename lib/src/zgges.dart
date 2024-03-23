import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zgeqrf.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zggbak.dart';
import 'package:lapack/src/zggbal.dart';
import 'package:lapack/src/zgghrd.dart';
import 'package:lapack/src/zhgeqz.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlascl.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/ztgsen.dart';
import 'package:lapack/src/zungqr.dart';
import 'package:lapack/src/zunmqr.dart';

void zgges(
  final String JOBVSL,
  final String JOBVSR,
  final String SORT,
  final bool Function(Complex, Complex) SELCTG,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> SDIM,
  final Array<Complex> ALPHA_,
  final Array<Complex> BETA_,
  final Matrix<Complex> VSL_,
  final int LDVSL,
  final Matrix<Complex> VSR_,
  final int LDVSR,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<bool> BWORK_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final VSL = VSL_.having(ld: LDVSL);
  final VSR = VSR_.having(ld: LDVSR);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final BWORK = BWORK_.having();
  final ALPHA = ALPHA_.having();
  final BETA = BETA_.having();

  const ZERO = 0.0, ONE = 1.0;
  bool CURSL, ILASCL, ILBSCL, ILVSL, ILVSR, LASTSL, LQUERY, WANTST;
  int I,
      ICOLS,
      IJOBVL,
      IJOBVR,
      ILEFT,
      IRIGHT,
      IROWS,
      IRWRK,
      ITAU,
      IWRK,
      LWKMIN,
      LWKOPT = 0;
  double ANRM, ANRMTO = 0, BIGNUM, BNRM, BNRMTO = 0, EPS, SMLNUM;
  final IDUM = Array<int>(1);
  final DIF = Array<double>(2);
  final IERR = Box(0), IHI = Box(0), ILO = Box(0);
  final PVSL = Box(0.0), PVSR = Box(0.0);

  // Decode the input arguments

  if (lsame(JOBVSL, 'N')) {
    IJOBVL = 1;
    ILVSL = false;
  } else if (lsame(JOBVSL, 'V')) {
    IJOBVL = 2;
    ILVSL = true;
  } else {
    IJOBVL = -1;
    ILVSL = false;
  }

  if (lsame(JOBVSR, 'N')) {
    IJOBVR = 1;
    ILVSR = false;
  } else if (lsame(JOBVSR, 'V')) {
    IJOBVR = 2;
    ILVSR = true;
  } else {
    IJOBVR = -1;
    ILVSR = false;
  }

  WANTST = lsame(SORT, 'S');

  // Test the input arguments

  INFO.value = 0;
  LQUERY = (LWORK == -1);
  if (IJOBVL <= 0) {
    INFO.value = -1;
  } else if (IJOBVR <= 0) {
    INFO.value = -2;
  } else if (!WANTST && !lsame(SORT, 'N')) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -5;
  } else if (LDA < max(1, N)) {
    INFO.value = -7;
  } else if (LDB < max(1, N)) {
    INFO.value = -9;
  } else if (LDVSL < 1 || (ILVSL && LDVSL < N)) {
    INFO.value = -14;
  } else if (LDVSR < 1 || (ILVSR && LDVSR < N)) {
    INFO.value = -16;
  }

  // Compute workspace
  //  (Note: Comments in the code beginning "Workspace:" describe the
  //   minimal amount of workspace needed at that point in the code,
  //   as well as the preferred amount for good performance.
  //   NB refers to the optimal block size for the immediately
  //   following subroutine, as returned by ILAENV.)

  if (INFO.value == 0) {
    LWKMIN = max(1, 2 * N);
    LWKOPT = max(1, N + N * ilaenv(1, 'ZGEQRF', ' ', N, 1, N, 0));
    LWKOPT = max(LWKOPT, N + N * ilaenv(1, 'ZUNMQR', ' ', N, 1, N, -1));
    if (ILVSL) {
      LWKOPT = max(LWKOPT, N + N * ilaenv(1, 'ZUNGQR', ' ', N, 1, N, -1));
    }
    WORK[1] = LWKOPT.toComplex();

    if (LWORK < LWKMIN && !LQUERY) INFO.value = -18;
  }

  if (INFO.value != 0) {
    xerbla('ZGGES ', -INFO.value);
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

  // Permute the matrix to make it more nearly triangular
  // (Real Workspace: need 6*N)

  ILEFT = 1;
  IRIGHT = N + 1;
  IRWRK = IRIGHT + N;
  zggbal('P', N, A, LDA, B, LDB, ILO, IHI, RWORK(ILEFT), RWORK(IRIGHT),
      RWORK(IRWRK), IERR);

  // Reduce B to triangular form (QR decomposition of B)
  // (Complex Workspace: need N, prefer N*NB)

  IROWS = IHI.value + 1 - ILO.value;
  ICOLS = N + 1 - ILO.value;
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

  // Initialize VSL
  // (Complex Workspace: need N, prefer N*NB)

  if (ILVSL) {
    zlaset('Full', N, N, Complex.zero, Complex.one, VSL, LDVSL);
    if (IROWS > 1) {
      zlacpy('L', IROWS - 1, IROWS - 1, B(ILO.value + 1, ILO.value), LDB,
          VSL(ILO.value + 1, ILO.value), LDVSL);
    }
    zungqr(IROWS, IROWS, IROWS, VSL(ILO.value, ILO.value), LDVSL, WORK(ITAU),
        WORK(IWRK), LWORK + 1 - IWRK, IERR);
  }

  // Initialize VSR

  if (ILVSR) zlaset('Full', N, N, Complex.zero, Complex.one, VSR, LDVSR);

  // Reduce to generalized Hessenberg form
  // (Workspace: none needed)

  zgghrd(JOBVSL, JOBVSR, N, ILO.value, IHI.value, A, LDA, B, LDB, VSL, LDVSL,
      VSR, LDVSR, IERR);

  SDIM.value = 0;

  // Perform QZ algorithm, computing Schur vectors if desired
  // (Complex Workspace: need N)
  // (Real Workspace: need N)

  IWRK = ITAU;
  zhgeqz(
      'S',
      JOBVSL,
      JOBVSR,
      N,
      ILO.value,
      IHI.value,
      A,
      LDA,
      B,
      LDB,
      ALPHA,
      BETA,
      VSL,
      LDVSL,
      VSR,
      LDVSR,
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
  } else {
    // Sort eigenvalues ALPHA/BETA if desired
    // (Workspace: none needed)

    if (WANTST) {
      // Undo scaling on eigenvalues before selecting

      if (ILASCL) {
        zlascl('G', 0, 0, ANRM, ANRMTO, N, 1, ALPHA.asMatrix(N), N, IERR);
      }
      if (ILBSCL) {
        zlascl('G', 0, 0, BNRM, BNRMTO, N, 1, BETA.asMatrix(N), N, IERR);
      }

      // Select eigenvalues

      for (I = 1; I <= N; I++) {
        BWORK[I] = SELCTG(ALPHA[I], BETA[I]);
      }

      ztgsen(
          0,
          ILVSL,
          ILVSR,
          BWORK,
          N,
          A,
          LDA,
          B,
          LDB,
          ALPHA,
          BETA,
          VSL,
          LDVSL,
          VSR,
          LDVSR,
          SDIM,
          PVSL,
          PVSR,
          DIF,
          WORK(IWRK),
          LWORK - IWRK + 1,
          IDUM,
          1,
          IERR);
      if (IERR.value == 1) INFO.value = N + 3;
    }

    // Apply back-permutation to VSL and VSR
    // (Workspace: none needed)

    if (ILVSL) {
      zggbak('P', 'L', N, ILO.value, IHI.value, RWORK(ILEFT), RWORK(IRIGHT), N,
          VSL, LDVSL, IERR);
    }
    if (ILVSR) {
      zggbak('P', 'R', N, ILO.value, IHI.value, RWORK(ILEFT), RWORK(IRIGHT), N,
          VSR, LDVSR, IERR);
    }

    // Undo scaling

    if (ILASCL) {
      zlascl('U', 0, 0, ANRMTO, ANRM, N, N, A, LDA, IERR);
      zlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA.asMatrix(N), N, IERR);
    }

    if (ILBSCL) {
      zlascl('U', 0, 0, BNRMTO, BNRM, N, N, B, LDB, IERR);
      zlascl('G', 0, 0, BNRMTO, BNRM, N, 1, BETA.asMatrix(N), N, IERR);
    }

    if (WANTST) {
      // Check if reordering is correct

      LASTSL = true;
      SDIM.value = 0;
      for (I = 1; I <= N; I++) {
        CURSL = SELCTG(ALPHA[I], BETA[I]);
        if (CURSL) SDIM.value++;
        if (CURSL && !LASTSL) INFO.value = N + 2;
        LASTSL = CURSL;
      }
    }
  }

  WORK[1] = LWKOPT.toComplex();
}
