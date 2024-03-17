import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
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
import 'package:lapack/src/dtgsen.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/dgeqrf.dart';
import 'package:lapack/src/xerbla.dart';

void dgges(
  final String JOBVSL,
  final String JOBVSR,
  final String SORT,
  final bool Function(double ZR, double ZI, double D) SELCTG,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Box<int> SDIM,
  final Array<double> ALPHAR_,
  final Array<double> ALPHAI_,
  final Array<double> BETA_,
  final Matrix<double> VSL_,
  final int LDVSL,
  final Matrix<double> VSR_,
  final int LDVSR,
  final Array<double> WORK_,
  final int LWORK,
  final Array<bool> BWORK_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final ALPHAR = ALPHAR_.having();
  final ALPHAI = ALPHAI_.having();
  final BETA = BETA_.having();
  final VSL = VSL_.having(ld: LDVSL);
  final VSR = VSR_.having(ld: LDVSR);
  final WORK = WORK_.having();
  final BWORK = BWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool CURSL, ILVSL, ILVSR, LASTSL, LQUERY, LST2SL, WANTST;
  int I,
      ICOLS,
      IJOBVL,
      IJOBVR,
      ILEFT,
      IP,
      IRIGHT,
      IROWS,
      ITAU,
      IWRK,
      MAXWRK = 0;
  double ANRMTO = 0, BNRMTO = 0;
  final IDUM = Array<int>(1);
  final DIF = Array<double>(2);
  final IHI = Box(0), ILO = Box(0), IERR = Box(0);
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
  } else if ((!WANTST) && (!lsame(SORT, 'N'))) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -5;
  } else if (LDA < max(1, N)) {
    INFO.value = -7;
  } else if (LDB < max(1, N)) {
    INFO.value = -9;
  } else if (LDVSL < 1 || (ILVSL && LDVSL < N)) {
    INFO.value = -15;
  } else if (LDVSR < 1 || (ILVSR && LDVSR < N)) {
    INFO.value = -17;
  }

  // Compute workspace
  // (Note: Comments in the code beginning "Workspace:" describe the
  // minimal amount of workspace needed at that point in the code,
  // as well as the preferred amount for good performance.
  // NB refers to the optimal block size for the immediately
  // following subroutine, as returned by ILAENV.)

  if (INFO.value == 0) {
    final int MINWRK;
    if (N > 0) {
      MINWRK = max(8 * N, 6 * N + 16);
      MAXWRK = MINWRK - N + N * ilaenv(1, 'DGEQRF', ' ', N, 1, N, 0);
      MAXWRK =
          max(MAXWRK, MINWRK - N + N * ilaenv(1, 'DORMQR', ' ', N, 1, N, -1));
      if (ILVSL) {
        MAXWRK =
            max(MAXWRK, MINWRK - N + N * ilaenv(1, 'DORGQR', ' ', N, 1, N, -1));
      }
    } else {
      MINWRK = 1;
      MAXWRK = 1;
    }
    WORK[1] = MAXWRK.toDouble();

    if (LWORK < MINWRK && !LQUERY) INFO.value = -19;
  }

  if (INFO.value != 0) {
    xerbla('DGGES ', -INFO.value);
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

  final EPS = dlamch('P');
  final SAFMIN = dlamch('S');
  final SAFMAX = ONE / SAFMIN;
  final SMLNUM = sqrt(SAFMIN) / EPS;
  final BIGNUM = ONE / SMLNUM;

  // Scale A if max element outside range [SMLNUM,BIGNUM]

  final ANRM = dlange('M', N, N, A, LDA, WORK);
  var ILASCL = false;
  if (ANRM > ZERO && ANRM < SMLNUM) {
    ANRMTO = SMLNUM;
    ILASCL = true;
  } else if (ANRM > BIGNUM) {
    ANRMTO = BIGNUM;
    ILASCL = true;
  }
  if (ILASCL) dlascl('G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR);

  // Scale B if max element outside range [SMLNUM,BIGNUM]

  final BNRM = dlange('M', N, N, B, LDB, WORK);
  var ILBSCL = false;
  if (BNRM > ZERO && BNRM < SMLNUM) {
    BNRMTO = SMLNUM;
    ILBSCL = true;
  } else if (BNRM > BIGNUM) {
    BNRMTO = BIGNUM;
    ILBSCL = true;
  }
  if (ILBSCL) dlascl('G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR);

  // Permute the matrix to make it more nearly triangular
  // (Workspace: need 6*N + 2*N space for storing balancing factors)

  ILEFT = 1;
  IRIGHT = N + 1;
  IWRK = IRIGHT + N;
  dggbal('P', N, A, LDA, B, LDB, ILO, IHI, WORK(ILEFT), WORK(IRIGHT),
      WORK(IWRK), IERR);

  // Reduce B to triangular form (QR decomposition of B)
  // (Workspace: need N, prefer N*NB)

  IROWS = IHI.value + 1 - ILO.value;
  ICOLS = N + 1 - ILO.value;
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

  // Initialize VSL
  // (Workspace: need N, prefer N*NB)

  if (ILVSL) {
    dlaset('Full', N, N, ZERO, ONE, VSL, LDVSL);
    if (IROWS > 1) {
      dlacpy('L', IROWS - 1, IROWS - 1, B(ILO.value + 1, ILO.value), LDB,
          VSL(ILO.value + 1, ILO.value), LDVSL);
    }
    dorgqr(IROWS, IROWS, IROWS, VSL(ILO.value, ILO.value), LDVSL, WORK(ITAU),
        WORK(IWRK), LWORK + 1 - IWRK, IERR);
  }

  // Initialize VSR

  if (ILVSR) dlaset('Full', N, N, ZERO, ONE, VSR, LDVSR);

  // Reduce to generalized Hessenberg form
  // (Workspace: none needed)

  dgghrd(JOBVSL, JOBVSR, N, ILO.value, IHI.value, A, LDA, B, LDB, VSL, LDVSL,
      VSR, LDVSR, IERR);

  // Perform QZ algorithm, computing Schur vectors if desired
  // (Workspace: need N)

  IWRK = ITAU;
  dhgeqz('S', JOBVSL, JOBVSR, N, ILO.value, IHI.value, A, LDA, B, LDB, ALPHAR,
      ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, WORK(IWRK), LWORK + 1 - IWRK, IERR);
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
    // (Workspace: need 4*N+16 )

    SDIM.value = 0;
    if (WANTST) {
      // Undo scaling on eigenvalues before SELCTGing

      if (ILASCL) {
        dlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAR.asMatrix(N), N, IERR);
        dlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAI.asMatrix(N), N, IERR);
      }
      if (ILBSCL) {
        dlascl('G', 0, 0, BNRMTO, BNRM, N, 1, BETA.asMatrix(N), N, IERR);
      }

      // Select eigenvalues

      for (I = 1; I <= N; I++) {
        BWORK[I] = SELCTG(ALPHAR[I], ALPHAI[I], BETA[I]);
      }

      dtgsen(
          0,
          ILVSL,
          ILVSR,
          BWORK,
          N,
          A,
          LDA,
          B,
          LDB,
          ALPHAR,
          ALPHAI,
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
      dggbak('P', 'L', N, ILO.value, IHI.value, WORK(ILEFT), WORK(IRIGHT), N,
          VSL, LDVSL, IERR);
    }

    if (ILVSR) {
      dggbak('P', 'R', N, ILO.value, IHI.value, WORK(ILEFT), WORK(IRIGHT), N,
          VSR, LDVSR, IERR);
    }

    // Check if unscaling would cause over/underflow, if so, rescale
    // (ALPHAR[I],ALPHAI[I],BETA[I]) so BETA[I] is on the order of
    // B[I][I] and ALPHAR[I] and ALPHAI[I] are on the order of A[I][I]

    if (ILASCL) {
      for (I = 1; I <= N; I++) {
        if (ALPHAI[I] != ZERO) {
          if ((ALPHAR[I] / SAFMAX) > (ANRMTO / ANRM) ||
              (SAFMIN / ALPHAR[I]) > (ANRM / ANRMTO)) {
            WORK[1] = (A[I][I] / ALPHAR[I]).abs();
            BETA[I] *= WORK[1];
            ALPHAR[I] *= WORK[1];
            ALPHAI[I] *= WORK[1];
          } else if ((ALPHAI[I] / SAFMAX) > (ANRMTO / ANRM) ||
              (SAFMIN / ALPHAI[I]) > (ANRM / ANRMTO)) {
            WORK[1] = (A[I][I + 1] / ALPHAI[I]).abs();
            BETA[I] *= WORK[1];
            ALPHAR[I] *= WORK[1];
            ALPHAI[I] *= WORK[1];
          }
        }
      }
    }

    if (ILBSCL) {
      for (I = 1; I <= N; I++) {
        if (ALPHAI[I] != ZERO) {
          if ((BETA[I] / SAFMAX) > (BNRMTO / BNRM) ||
              (SAFMIN / BETA[I]) > (BNRM / BNRMTO)) {
            WORK[1] = (B[I][I] / BETA[I]).abs();
            BETA[I] *= WORK[1];
            ALPHAR[I] *= WORK[1];
            ALPHAI[I] *= WORK[1];
          }
        }
      }
    }

    // Undo scaling

    if (ILASCL) {
      dlascl('H', 0, 0, ANRMTO, ANRM, N, N, A, LDA, IERR);
      dlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAR.asMatrix(N), N, IERR);
      dlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAI.asMatrix(N), N, IERR);
    }

    if (ILBSCL) {
      dlascl('U', 0, 0, BNRMTO, BNRM, N, N, B, LDB, IERR);
      dlascl('G', 0, 0, BNRMTO, BNRM, N, 1, BETA.asMatrix(N), N, IERR);
    }

    if (WANTST) {
      // Check if reordering is correct

      LASTSL = true;
      LST2SL = true;
      SDIM.value = 0;
      IP = 0;
      for (I = 1; I <= N; I++) {
        CURSL = SELCTG(ALPHAR[I], ALPHAI[I], BETA[I]);
        if (ALPHAI[I] == ZERO) {
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
  }

  WORK[1] = MAXWRK.toDouble();
}
