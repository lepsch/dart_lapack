import 'dart:math';

import 'package:lapack/src/blas/drot.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dhgeqz.dart';
import 'package:lapack/src/dlanhs.dart';
import 'package:lapack/src/dlaqz3.dart';
import 'package:lapack/src/dlaqz4.dart';
import 'package:lapack/src/dlartg.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlaqz0(
  final String WANTS,
  final String WANTQ,
  final String WANTZ,
  final int N,
  final int ILO,
  final int IHI,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Array<double> ALPHAR_,
  final Array<double> ALPHAI_,
  final Array<double> BETA_,
  final Matrix<double> Q_,
  final int LDQ,
  final Matrix<double> Z_,
  final int LDZ,
  final Array<double> WORK_,
  final int LWORK,
  final int REC,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final ALPHAR = ALPHAR_.having();
  final ALPHAI = ALPHAI_.having();
  final BETA = BETA_.having();
  final Q = Q_.having(ld: LDQ);
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0; //, HALF = 0.5;

  double SMLNUM,
      ULP,
      ESHIFT = 0,
      SAFMIN,
      // SAFMAX,
      SWAP,
      BNORM,
      BTOL;
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
      RCOST,
      I;
  bool ILSCHUR = false, ILQ = false, ILZ = false;
  String JBCMPZ;
  final AED_INFO = Box(0),
      SWEEP_INFO = Box(0),
      NORM_INFO = Box(0),
      N_UNDEFLATED = Box(0),
      N_DEFLATED = Box(0);
  final C1 = Box(0.0), S1 = Box(0.0), TEMP = Box(0.0);

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
    xerbla('DLAQZ0', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N <= 0) {
    WORK[1] = 1.toDouble();
    return;
  }

  // Get the parameters

  JBCMPZ = '${WANTS[0]}${WANTQ[0]}${WANTZ[0]}';

  NMIN = ilaenv(12, 'DLAQZ0', JBCMPZ, N, ILO, IHI, LWORK);

  NWR = ilaenv(13, 'DLAQZ0', JBCMPZ, N, ILO, IHI, LWORK);
  NWR = max(2, NWR);
  NWR = min(IHI - ILO + 1, min((N - 1) ~/ 3, NWR));

  NIBBLE = ilaenv(14, 'DLAQZ0', JBCMPZ, N, ILO, IHI, LWORK);

  NSR = ilaenv(15, 'DLAQZ0', JBCMPZ, N, ILO, IHI, LWORK);
  NSR = min(NSR, min((N + 6) ~/ 9, IHI - ILO));
  NSR = max(2, NSR - (NSR % 2));

  RCOST = ilaenv(17, 'DLAQZ0', JBCMPZ, N, ILO, IHI, LWORK);
  ITEMP1 = NSR ~/ sqrt(1 + 2 * NSR / (RCOST.toDouble() / 100 * N));
  ITEMP1 = ((ITEMP1 - 1) ~/ 4) * 4 + 4;
  NBR = NSR + ITEMP1;

  if (N < NMIN || REC >= 2) {
    dhgeqz(WANTS, WANTQ, WANTZ, N, ILO, IHI, A, LDA, B, LDB, ALPHAR, ALPHAI,
        BETA, Q, LDQ, Z, LDZ, WORK, LWORK, INFO);
    return;
  }

  // Find out required workspace

  // Workspace query to dlaqz3
  NW = max(NWR, NMIN);
  dlaqz3(
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
      ALPHAR,
      ALPHAI,
      BETA,
      WORK.asMatrix(NW),
      NW,
      WORK.asMatrix(NW),
      NW,
      WORK,
      -1,
      REC,
      AED_INFO);
  ITEMP1 = WORK[1].toInt();
  // Workspace query to dlaqz4
  dlaqz4(
      ILSCHUR,
      ILQ,
      ILZ,
      N,
      ILO,
      IHI,
      NSR,
      NBR,
      ALPHAR,
      ALPHAI,
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
    WORK[1] = LWORKREQ.toDouble();
    return;
  } else if (LWORK < LWORKREQ) {
    INFO.value = -19;
  }
  if (INFO.value != 0) {
    xerbla('DLAQZ0', INFO.value);
    return;
  }

  // Initialize Q and Z

  if (IWANTQ == 3) dlaset('FULL', N, N, ZERO, ONE, Q, LDQ);
  if (IWANTZ == 3) dlaset('FULL', N, N, ZERO, ONE, Z, LDZ);

  // Get machine constants
  SAFMIN = dlamch('SAFE MINIMUM');
  // SAFMAX = ONE / SAFMIN;
  ULP = dlamch('PRECISION');
  SMLNUM = SAFMIN * (N.toDouble() / ULP);

  BNORM = dlanhs('F', IHI - ILO + 1, B(ILO, ILO), LDB, WORK);
  BTOL = max(SAFMIN, ULP * BNORM);

  ISTART = ILO;
  ISTOP = IHI;
  MAXIT = 3 * (IHI - ILO + 1);
  LD = 0;

  for (IITER = 1; IITER <= MAXIT; IITER++) {
    if (IITER >= MAXIT) {
      INFO.value = ISTOP + 1;
      break;
    }
    if (ISTART + 1 >= ISTOP) {
      ISTOP = ISTART;
      break;
    }

    // Check deflations at the end
    if ((A[ISTOP - 1][ISTOP - 2]).abs() <=
        max(
          SMLNUM,
          ULP *
              ((A[ISTOP - 1][ISTOP - 1]).abs() +
                  (A[ISTOP - 2][ISTOP - 2]).abs()),
        )) {
      A[ISTOP - 1][ISTOP - 2] = ZERO;
      ISTOP = ISTOP - 2;
      LD = 0;
      ESHIFT = ZERO;
    } else if ((A[ISTOP][ISTOP - 1]).abs() <=
        max(
          SMLNUM,
          ULP * ((A[ISTOP][ISTOP]).abs() + (A[ISTOP - 1][ISTOP - 1]).abs()),
        )) {
      A[ISTOP][ISTOP - 1] = ZERO;
      ISTOP = ISTOP - 1;
      LD = 0;
      ESHIFT = ZERO;
    }
    // Check deflations at the start
    if ((A[ISTART + 2][ISTART + 1]).abs() <=
        max(
          SMLNUM,
          ULP *
              ((A[ISTART + 1][ISTART + 1]).abs() +
                  (A[ISTART + 2][ISTART + 2]).abs()),
        )) {
      A[ISTART + 2][ISTART + 1] = ZERO;
      ISTART = ISTART + 2;
      LD = 0;
      ESHIFT = ZERO;
    } else if ((A[ISTART + 1][ISTART]).abs() <=
        max(
          SMLNUM,
          ULP * ((A[ISTART][ISTART]).abs() + (A[ISTART + 1][ISTART + 1]).abs()),
        )) {
      A[ISTART + 1][ISTART] = ZERO;
      ISTART = ISTART + 1;
      LD = 0;
      ESHIFT = ZERO;
    }

    if (ISTART + 1 >= ISTOP) {
      break;
    }

    // Check interior deflations
    ISTART2 = ISTART;
    for (K = ISTOP; K >= ISTART + 1; K--) {
      if ((A[K][K - 1]).abs() <=
          max(SMLNUM, ULP * ((A[K][K]).abs() + (A[K - 1][K - 1]).abs()))) {
        A[K][K - 1] = ZERO;
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
          dlartg(B[K2 - 1][K2], B[K2 - 1][K2 - 1], C1, S1, TEMP);
          B[K2 - 1][K2] = TEMP.value;
          B[K2 - 1][K2 - 1] = ZERO;
          drot(K2 - 2 - ISTARTM + 1, B(ISTARTM, K2).asArray(), 1,
              B(ISTARTM, K2 - 1).asArray(), 1, C1.value, S1.value);
          drot(min(K2 + 1, ISTOP) - ISTARTM + 1, A(ISTARTM, K2).asArray(), 1,
              A(ISTARTM, K2 - 1).asArray(), 1, C1.value, S1.value);
          if (ILZ) {
            drot(N, Z(1, K2).asArray(), 1, Z(1, K2 - 1).asArray(), 1, C1.value,
                S1.value);
          }

          if (K2 < ISTOP) {
            dlartg(A[K2][K2 - 1], A[K2 + 1][K2 - 1], C1, S1, TEMP);
            A[K2][K2 - 1] = TEMP.value;
            A[K2 + 1][K2 - 1] = ZERO;
            drot(ISTOPM - K2 + 1, A(K2, K2).asArray(), LDA,
                A(K2 + 1, K2).asArray(), LDA, C1.value, S1.value);
            drot(ISTOPM - K2 + 1, B(K2, K2).asArray(), LDB,
                B(K2 + 1, K2).asArray(), LDB, C1.value, S1.value);
            if (ILQ) {
              drot(N, Q(1, K2).asArray(), 1, Q(1, K2 + 1).asArray(), 1,
                  C1.value, S1.value);
            }
          }
        }

        if (ISTART2 < ISTOP) {
          dlartg(A[ISTART2][ISTART2], A[ISTART2 + 1][ISTART2], C1, S1, TEMP);
          A[ISTART2][ISTART2] = TEMP.value;
          A[ISTART2 + 1][ISTART2] = ZERO;
          drot(
              ISTOPM - (ISTART2 + 1) + 1,
              A(ISTART2, ISTART2 + 1).asArray(),
              LDA,
              A(ISTART2 + 1, ISTART2 + 1).asArray(),
              LDA,
              C1.value,
              S1.value);
          drot(
              ISTOPM - (ISTART2 + 1) + 1,
              B(ISTART2, ISTART2 + 1).asArray(),
              LDB,
              B(ISTART2 + 1, ISTART2 + 1).asArray(),
              LDB,
              C1.value,
              S1.value);
          if (ILQ) {
            drot(N, Q(1, ISTART2).asArray(), 1, Q(1, ISTART2 + 1).asArray(), 1,
                C1.value, S1.value);
          }
        }

        ISTART2 = ISTART2 + 1;
      }
      K = K - 1;
    }

    // istart2 now points to the top of the bottom right
    // unreduced Hessenberg block
    if (ISTART2 >= ISTOP) {
      ISTOP = ISTART2 - 1;
      LD = 0;
      ESHIFT = ZERO;
      continue;
    }

    NW = NWR;
    NSHIFTS = NSR;
    NBLOCK = NBR;

    if (ISTOP - ISTART2 + 1 < NMIN) {
      // Setting nw to the size of the subblock will make AED deflate
      // all the eigenvalues. This is slightly more efficient than just
      // using DHGEQZ because the off diagonal part gets updated via BLAS.
      if (ISTOP - ISTART + 1 < NMIN) {
        NW = ISTOP - ISTART + 1;
        ISTART2 = ISTART;
      } else {
        NW = ISTOP - ISTART2 + 1;
      }
    }

    // Time for AED

    dlaqz3(
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
        ALPHAR,
        ALPHAI,
        BETA,
        WORK.asMatrix(NW),
        NW,
        WORK(pow(NW, 2).toInt() + 1).asMatrix(NW),
        NW,
        WORK(2 * pow(NW, 2).toInt() + 1),
        LWORK - 2 * pow(NW, 2).toInt(),
        REC,
        AED_INFO);

    if (N_DEFLATED.value > 0) {
      ISTOP = ISTOP - N_DEFLATED.value;
      LD = 0;
      ESHIFT = ZERO;
    }
    if (100 * N_DEFLATED.value >
            NIBBLE * (N_DEFLATED.value + N_UNDEFLATED.value) ||
        ISTOP - ISTART2 + 1 < NMIN) {
      // AED has uncovered many eigenvalues. Skip a QZ sweep and run
      // AED again.
      continue;
    }

    LD = LD + 1;

    NS = min(NSHIFTS, ISTOP - ISTART2);
    NS = min(NS, N_UNDEFLATED.value);
    SHIFTPOS = ISTOP - N_UNDEFLATED.value + 1;

    // Shuffle shifts to put double shifts in front
    // This ensures that we don't split up a double shift

    for (I = SHIFTPOS;
        2 < 0
            ? I >= SHIFTPOS + N_UNDEFLATED.value - 1
            : I <= SHIFTPOS + N_UNDEFLATED.value - 1;
        I += 2) {
      if (ALPHAI[I] != -ALPHAI[I + 1]) {
        SWAP = ALPHAR[I];
        ALPHAR[I] = ALPHAR[I + 1];
        ALPHAR[I + 1] = ALPHAR[I + 2];
        ALPHAR[I + 2] = SWAP;

        SWAP = ALPHAI[I];
        ALPHAI[I] = ALPHAI[I + 1];
        ALPHAI[I + 1] = ALPHAI[I + 2];
        ALPHAI[I + 2] = SWAP;

        SWAP = BETA[I];
        BETA[I] = BETA[I + 1];
        BETA[I + 1] = BETA[I + 2];
        BETA[I + 2] = SWAP;
      }
    }

    if ((LD % 6) == 0) {
      // Exceptional shift.  Chosen for no particularly good reason.

      if ((MAXIT.toDouble() * SAFMIN) * (A[ISTOP][ISTOP - 1]).abs() <
          (A[ISTOP - 1][ISTOP - 1]).abs()) {
        ESHIFT = A[ISTOP][ISTOP - 1] / B[ISTOP - 1][ISTOP - 1];
      } else {
        ESHIFT = ESHIFT + ONE / (SAFMIN * MAXIT.toDouble());
      }
      ALPHAR[SHIFTPOS] = ONE;
      ALPHAR[SHIFTPOS + 1] = ZERO;
      ALPHAI[SHIFTPOS] = ZERO;
      ALPHAI[SHIFTPOS + 1] = ZERO;
      BETA[SHIFTPOS] = ESHIFT;
      BETA[SHIFTPOS + 1] = ESHIFT;
      NS = 2;
    }

    // Time for a QZ sweep

    dlaqz4(
        ILSCHUR,
        ILQ,
        ILZ,
        N,
        ISTART2,
        ISTOP,
        NS,
        NBLOCK,
        ALPHAR(SHIFTPOS),
        ALPHAI(SHIFTPOS),
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

  // Call DHGEQZ to normalize the eigenvalue blocks and set the eigenvalues
  // If all the eigenvalues have been found, DHGEQZ will not do any iterations
  // and only normalize the blocks. In case of a rare convergence failure,
  // the single shift might perform better.

  dhgeqz(WANTS, WANTQ, WANTZ, N, ILO, IHI, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA,
      Q, LDQ, Z, LDZ, WORK, LWORK, NORM_INFO);

  INFO.value = NORM_INFO.value;
}
