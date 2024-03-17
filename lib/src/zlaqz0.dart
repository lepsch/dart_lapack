import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhgeqz.dart';
import 'package:lapack/src/zlanhs.dart';
import 'package:lapack/src/zlaqz2.dart';
import 'package:lapack/src/zlaqz3.dart';
import 'package:lapack/src/zlartg.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zrot.dart';

void zlaqz0(
  final String WANTS,
  final String WANTQ,
  final String WANTZ,
  final int N,
  final int ILO,
  final int IHI,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<Complex> ALPHA_,
  final Array<Complex> BETA_,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final int REC,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final ALPHA = ALPHA_.having();
  final BETA = BETA_.having();
  final Q = Q_.having(ld: LDQ);
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  // const ZERO = 0.0, ONE = 1.0, HALF = 0.5;
  double SMLNUM,
      ULP,
      SAFMIN,
      // SAFMAX,
      // TEMPR,
      BNORM,
      BTOL;
  Complex ESHIFT = Complex.zero;
  int ISTART,
      ISTOP,
      IITER,
      MAXIT,
      ISTART2,
      K,
      LD,
      NSHIFTS,
      NBLOCK,
      NW,
      NMIN,
      NIBBLE,
      NS,
      SHIFTPOS,
      LWORKREQ,
      K2,
      ISTARTM,
      ISTOPM,
      IWANTS,
      IWANTQ,
      IWANTZ,
      NWR,
      NBR,
      NSR,
      ITEMP1,
      ITEMP2,
      RCOST;
  bool ILSCHUR = false, ILQ = false, ILZ = false;
  String JBCMPZ = '';
  final N_UNDEFLATED = Box(0),
      N_DEFLATED = Box(0),
      AED_INFO = Box(0),
      NORM_INFO = Box(0),
      SWEEP_INFO = Box(0);
  final C1 = Box(0.0);
  final S1 = Box(Complex.zero), TEMP = Box(Complex.zero);

  // Decode wantS,wantQ,wantZ

  if (lsame(WANTS, 'E')) {
    ILSCHUR = false;
    IWANTS = 1;
  } else if (lsame(WANTS, 'S')) {
    ILSCHUR = true;
    IWANTS = 2;
  } else {
    IWANTS = 0;
  }

  if (lsame(WANTQ, 'N')) {
    ILQ = false;
    IWANTQ = 1;
  } else if (lsame(WANTQ, 'V')) {
    ILQ = true;
    IWANTQ = 2;
  } else if (lsame(WANTQ, 'I')) {
    ILQ = true;
    IWANTQ = 3;
  } else {
    IWANTQ = 0;
  }

  if (lsame(WANTZ, 'N')) {
    ILZ = false;
    IWANTZ = 1;
  } else if (lsame(WANTZ, 'V')) {
    ILZ = true;
    IWANTZ = 2;
  } else if (lsame(WANTZ, 'I')) {
    ILZ = true;
    IWANTZ = 3;
  } else {
    IWANTZ = 0;
  }

  // Check Argument Values

  INFO.value = 0;
  if (IWANTS == 0) {
    INFO.value = -1;
  } else if (IWANTQ == 0) {
    INFO.value = -2;
  } else if (IWANTZ == 0) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (ILO < 1) {
    INFO.value = -5;
  } else if (IHI > N || IHI < ILO - 1) {
    INFO.value = -6;
  } else if (LDA < N) {
    INFO.value = -8;
  } else if (LDB < N) {
    INFO.value = -10;
  } else if (LDQ < 1 || (ILQ && LDQ < N)) {
    INFO.value = -15;
  } else if (LDZ < 1 || (ILZ && LDZ < N)) {
    INFO.value = -17;
  }
  if (INFO.value != 0) {
    xerbla('ZLAQZ0', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N <= 0) {
    WORK[1] = 1.toComplex();
    return;
  }

  // Get the parameters

  JBCMPZ = '$WANTS$WANTQ$WANTZ';

  NMIN = ilaenv(12, 'ZLAQZ0', JBCMPZ, N, ILO, IHI, LWORK);

  NWR = ilaenv(13, 'ZLAQZ0', JBCMPZ, N, ILO, IHI, LWORK);
  NWR = max(2, NWR);
  NWR = min(IHI - ILO + 1, min((N - 1) ~/ 3, NWR));

  NIBBLE = ilaenv(14, 'ZLAQZ0', JBCMPZ, N, ILO, IHI, LWORK);

  NSR = ilaenv(15, 'ZLAQZ0', JBCMPZ, N, ILO, IHI, LWORK);
  NSR = min(NSR, min((N + 6) ~/ 9, IHI - ILO));
  NSR = max(2, NSR - (NSR % 2));

  RCOST = ilaenv(17, 'ZLAQZ0', JBCMPZ, N, ILO, IHI, LWORK);
  ITEMP1 = NSR ~/ sqrt(1 + 2 * NSR / (RCOST.toDouble() / 100 * N));
  ITEMP1 = ((ITEMP1 - 1) ~/ 4) * 4 + 4;
  NBR = NSR + ITEMP1;

  if (N < NMIN || REC >= 2) {
    zhgeqz(WANTS, WANTQ, WANTZ, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, Q,
        LDQ, Z, LDZ, WORK, LWORK, RWORK, INFO);
    return;
  }

  // Find out required workspace

  // Workspace query to ZLAQZ2
  NW = max(NWR, NMIN);
  zlaqz2(
      ILSCHUR,
      ILQ,
      ILZ,
      N,
      ILO,
      IHI,
      NW,
      A,
      LDA,
      B,
      LDB,
      Q,
      LDQ,
      Z,
      LDZ,
      N_UNDEFLATED,
      N_DEFLATED,
      ALPHA,
      BETA,
      WORK.asMatrix(NW),
      NW,
      WORK.asMatrix(NW),
      NW,
      WORK,
      -1,
      RWORK,
      REC,
      AED_INFO);
  ITEMP1 = WORK[1].toInt();
  // Workspace query to ZLAQZ3
  zlaqz3(
      ILSCHUR,
      ILQ,
      ILZ,
      N,
      ILO,
      IHI,
      NSR,
      NBR,
      ALPHA,
      BETA,
      A,
      LDA,
      B,
      LDB,
      Q,
      LDQ,
      Z,
      LDZ,
      WORK.asMatrix(NBR),
      NBR,
      WORK.asMatrix(NBR),
      NBR,
      WORK,
      -1,
      SWEEP_INFO);
  ITEMP2 = WORK[1].toInt();

  LWORKREQ =
      max(ITEMP1 + 2 * pow(NW, 2).toInt(), ITEMP2 + 2 * pow(NBR, 2).toInt());
  if (LWORK == -1) {
    WORK[1] = LWORKREQ.toComplex();
    return;
  } else if (LWORK < LWORKREQ) {
    INFO.value = -19;
  }
  if (INFO.value != 0) {
    xerbla('ZLAQZ0', INFO.value);
    return;
  }

  // Initialize Q and Z

  if (IWANTQ == 3) zlaset('FULL', N, N, Complex.zero, Complex.one, Q, LDQ);
  if (IWANTZ == 3) zlaset('FULL', N, N, Complex.zero, Complex.one, Z, LDZ);

  // Get machine constants
  SAFMIN = dlamch('SAFE MINIMUM');
  // SAFMAX = ONE/SAFMIN;
  ULP = dlamch('PRECISION');
  SMLNUM = SAFMIN * (N.toDouble() / ULP);

  BNORM = zlanhs('F', IHI - ILO + 1, B(ILO, ILO), LDB, RWORK);
  BTOL = max(SAFMIN, ULP * BNORM);

  ISTART = ILO;
  ISTOP = IHI;
  MAXIT = 30 * (IHI - ILO + 1);
  LD = 0;

  for (IITER = 1; IITER <= MAXIT; IITER++) {
    if (IITER >= MAXIT) {
      INFO.value = ISTOP + 1;
      INFO.value = NORM_INFO.value;
      return;
    }
    if (ISTART + 1 >= ISTOP) {
      ISTOP = ISTART;
      break;
    }

    // Check deflations at the end
    if (A[ISTOP][ISTOP - 1].abs() <=
        max(
            SMLNUM,
            ULP *
                ((A[ISTOP][ISTOP]).abs() + (A[ISTOP - 1][ISTOP - 1]).abs()))) {
      A[ISTOP][ISTOP - 1] = Complex.zero;
      ISTOP--;
      LD = 0;
      ESHIFT = Complex.zero;
    }
    // Check deflations at the start
    if ((A[ISTART + 1][ISTART]).abs() <=
        max(
            SMLNUM,
            ULP *
                ((A[ISTART][ISTART]).abs() +
                    (A[ISTART + 1][ISTART + 1]).abs()))) {
      A[ISTART + 1][ISTART] = Complex.zero;
      ISTART++;
      LD = 0;
      ESHIFT = Complex.zero;
    }

    if (ISTART + 1 >= ISTOP) {
      break;
    }

    // Check interior deflations
    ISTART2 = ISTART;
    for (K = ISTOP; K >= ISTART + 1; K--) {
      if (A[K][K - 1].abs() <=
          max(SMLNUM, ULP * ((A[K][K]).abs() + (A[K - 1][K - 1]).abs()))) {
        A[K][K - 1] = Complex.zero;
        ISTART2 = K;
        break;
      }
    }

    // Get range to apply rotations to
    if (ILSCHUR) {
      ISTARTM = 1;
      ISTOPM = N;
    } else {
      ISTARTM = ISTART2;
      ISTOPM = ISTOP;
    }

    // Check infinite eigenvalues, this is done without blocking so might
    // slow down the method when many infinite eigenvalues are present
    K = ISTOP;
    while (K >= ISTART2) {
      if ((B[K][K]).abs() < BTOL) {
        // A diagonal element of B is negligible, move it
        // to the top and deflate it

        for (K2 = K; K2 >= ISTART2 + 1; K2--) {
          zlartg(B[K2 - 1][K2], B[K2 - 1][K2 - 1], C1, S1, TEMP);
          B[K2 - 1][K2] = TEMP.value;
          B[K2 - 1][K2 - 1] = Complex.zero;
          zrot(K2 - 2 - ISTARTM + 1, B(ISTARTM, K2).asArray(), 1,
              B(ISTARTM, K2 - 1).asArray(), 1, C1.value, S1.value);
          zrot(min(K2 + 1, ISTOP) - ISTARTM + 1, A(ISTARTM, K2).asArray(), 1,
              A(ISTARTM, K2 - 1).asArray(), 1, C1.value, S1.value);
          if (ILZ) {
            zrot(N, Z(1, K2).asArray(), 1, Z(1, K2 - 1).asArray(), 1, C1.value,
                S1.value);
          }

          if (K2 < ISTOP) {
            zlartg(A[K2][K2 - 1], A[K2 + 1][K2 - 1], C1, S1, TEMP);
            A[K2][K2 - 1] = TEMP.value;
            A[K2 + 1][K2 - 1] = Complex.zero;
            zrot(ISTOPM - K2 + 1, A(K2, K2).asArray(), LDA,
                A(K2 + 1, K2).asArray(), LDA, C1.value, S1.value);
            zrot(ISTOPM - K2 + 1, B(K2, K2).asArray(), LDB,
                B(K2 + 1, K2).asArray(), LDB, C1.value, S1.value);
            if (ILQ) {
              zrot(N, Q(1, K2).asArray(), 1, Q(1, K2 + 1).asArray(), 1,
                  C1.value, S1.value.conjugate());
            }
          }
        }

        if (ISTART2 < ISTOP) {
          zlartg(A[ISTART2][ISTART2], A[ISTART2 + 1][ISTART2], C1, S1, TEMP);
          A[ISTART2][ISTART2] = TEMP.value;
          A[ISTART2 + 1][ISTART2] = Complex.zero;
          zrot(
              ISTOPM - (ISTART2 + 1) + 1,
              A(ISTART2, ISTART2 + 1).asArray(),
              LDA,
              A(ISTART2 + 1, ISTART2 + 1).asArray(),
              LDA,
              C1.value,
              S1.value);
          zrot(
              ISTOPM - (ISTART2 + 1) + 1,
              B(ISTART2, ISTART2 + 1).asArray(),
              LDB,
              B(ISTART2 + 1, ISTART2 + 1).asArray(),
              LDB,
              C1.value,
              S1.value);
          if (ILQ) {
            zrot(N, Q(1, ISTART2).asArray(), 1, Q(1, ISTART2 + 1).asArray(), 1,
                C1.value, S1.value.conjugate());
          }
        }

        ISTART2++;
      }
      K--;
    }

    // istart2 now points to the top of the bottom right
    // unreduced Hessenberg block
    if (ISTART2 >= ISTOP) {
      ISTOP = ISTART2 - 1;
      LD = 0;
      ESHIFT = Complex.zero;
      continue;
    }

    NW = NWR;
    NSHIFTS = NSR;
    NBLOCK = NBR;

    if (ISTOP - ISTART2 + 1 < NMIN) {
      // Setting nw to the size of the subblock will make AED deflate
      // all the eigenvalues. This is slightly more efficient than just
      // using qz_small because the off diagonal part gets updated via BLAS.
      if (ISTOP - ISTART + 1 < NMIN) {
        NW = ISTOP - ISTART + 1;
        ISTART2 = ISTART;
      } else {
        NW = ISTOP - ISTART2 + 1;
      }
    }

    // Time for AED

    zlaqz2(
        ILSCHUR,
        ILQ,
        ILZ,
        N,
        ISTART2,
        ISTOP,
        NW,
        A,
        LDA,
        B,
        LDB,
        Q,
        LDQ,
        Z,
        LDZ,
        N_UNDEFLATED,
        N_DEFLATED,
        ALPHA,
        BETA,
        WORK.asMatrix(NW),
        NW,
        WORK(pow(NW, 2).toInt() + 1).asMatrix(NW),
        NW,
        WORK(2 * pow(NW, 2).toInt() + 1),
        LWORK - 2 * pow(NW, 2).toInt(),
        RWORK,
        REC,
        AED_INFO);

    if (N_DEFLATED.value > 0) {
      ISTOP -= N_DEFLATED.value;
      LD = 0;
      ESHIFT = Complex.zero;
    }
    if (100 * N_DEFLATED.value >
            NIBBLE * (N_DEFLATED.value + N_UNDEFLATED.value) ||
        ISTOP - ISTART2 + 1 < NMIN) {
      // AED has uncovered many eigenvalues. Skip a QZ sweep and run
      // AED again.
      continue;
    }

    LD++;

    NS = min(NSHIFTS, ISTOP - ISTART2);
    NS = min(NS, N_UNDEFLATED.value);
    SHIFTPOS = ISTOP - N_UNDEFLATED.value + 1;

    if ((LD % 6) == 0) {
      // Exceptional shift.  Chosen for no particularly good reason.

      if ((MAXIT.toDouble() * SAFMIN) * (A[ISTOP][ISTOP - 1]).abs() <
          (A[ISTOP - 1][ISTOP - 1]).abs()) {
        ESHIFT = A[ISTOP][ISTOP - 1] / B[ISTOP - 1][ISTOP - 1];
      } else {
        ESHIFT += Complex.one / (SAFMIN * MAXIT).toComplex();
      }
      ALPHA[SHIFTPOS] = Complex.one;
      BETA[SHIFTPOS] = ESHIFT;
      NS = 1;
    }

    // Time for a QZ sweep

    zlaqz3(
        ILSCHUR,
        ILQ,
        ILZ,
        N,
        ISTART2,
        ISTOP,
        NS,
        NBLOCK,
        ALPHA(SHIFTPOS),
        BETA(SHIFTPOS),
        A,
        LDA,
        B,
        LDB,
        Q,
        LDQ,
        Z,
        LDZ,
        WORK.asMatrix(NBLOCK),
        NBLOCK,
        WORK(pow(NBLOCK, 2).toInt() + 1).asMatrix(NBLOCK),
        NBLOCK,
        WORK(2 * pow(NBLOCK, 2).toInt() + 1),
        LWORK - 2 * pow(NBLOCK, 2).toInt(),
        SWEEP_INFO);
  }

  // Call ZHGEQZ to normalize the eigenvalue blocks and set the eigenvalues
  // If all the eigenvalues have been found, ZHGEQZ will not do any iterations
  // and only normalize the blocks. In case of a rare convergence failure,
  // the single shift might perform better.

  zhgeqz(WANTS, WANTQ, WANTZ, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, Q, LDQ,
      Z, LDZ, WORK, LWORK, RWORK, NORM_INFO);

  INFO.value = NORM_INFO.value;
}
