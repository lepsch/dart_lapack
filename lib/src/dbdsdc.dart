import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlanst.dart';
import 'package:lapack/src/dlartg.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dlasd0.dart';
import 'package:lapack/src/dlasda.dart';
import 'package:lapack/src/dlasdq.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dlasr.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dbdsdc(
  final String UPLO,
  final String COMPQ,
  final int N,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<double> U_,
  final int LDU,
  final Matrix<double> VT_,
  final int LDVT,
  final Array<double> Q_,
  final Array<int> IQ_,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final E = E_.having();
  final U = U_.having(ld: LDU);
  final VT = VT_.having(ld: LDVT);
  final Q = Q_.having();
  final IQ = IQ_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  int DIFL = 0,
      DIFR = 0,
      GIVCOL = 0,
      GIVNUM = 0,
      GIVPTR = 0,
      I,
      IC = 0,
      ICOMPQ,
      II,
      IS = 0,
      IU = 0,
      IUPLO,
      IVT = 0,
      J,
      K = 0,
      KK,
      MLVL,
      NM1,
      NSIZE,
      PERM = 0,
      POLES = 0,
      QSTART,
      SMLSIZ,
      SMLSZP,
      SQRE,
      START,
      WSTART,
      Z = 0;
  double EPS, ORGNRM, P;
  final CS = Box(0.0), R = Box(0.0), SN = Box(0.0);
  final IERR = Box(0);

  // Test the input parameters.

  INFO.value = 0;

  IUPLO = 0;
  if (lsame(UPLO, 'U')) IUPLO = 1;
  if (lsame(UPLO, 'L')) IUPLO = 2;
  if (lsame(COMPQ, 'N')) {
    ICOMPQ = 0;
  } else if (lsame(COMPQ, 'P')) {
    ICOMPQ = 1;
  } else if (lsame(COMPQ, 'I')) {
    ICOMPQ = 2;
  } else {
    ICOMPQ = -1;
  }
  if (IUPLO == 0) {
    INFO.value = -1;
  } else if (ICOMPQ < 0) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if ((LDU < 1) || ((ICOMPQ == 2) && (LDU < N))) {
    INFO.value = -7;
  } else if ((LDVT < 1) || ((ICOMPQ == 2) && (LDVT < N))) {
    INFO.value = -9;
  }
  if (INFO.value != 0) {
    xerbla('DBDSDC', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;
  SMLSIZ = ilaenv(9, 'DBDSDC', ' ', 0, 0, 0, 0);
  if (N == 1) {
    if (ICOMPQ == 1) {
      Q[1] = sign(ONE, D[1]).toDouble();
      Q[1 + SMLSIZ * N] = ONE;
    } else if (ICOMPQ == 2) {
      U[1][1] = sign(ONE, D[1]).toDouble();
      VT[1][1] = ONE;
    }
    D[1] = (D[1]).abs();
    return;
  }
  NM1 = N - 1;

  // If matrix lower bidiagonal, rotate to be upper bidiagonal
  // by applying Givens rotations on the left

  WSTART = 1;
  QSTART = 3;
  if (ICOMPQ == 1) {
    dcopy(N, D, 1, Q(1), 1);
    dcopy(N - 1, E, 1, Q(N + 1), 1);
  }
  if (IUPLO == 2) {
    QSTART = 5;
    if (ICOMPQ == 2) WSTART = 2 * N - 1;
    for (I = 1; I <= N - 1; I++) {
      dlartg(D[I], E[I], CS, SN, R);
      D[I] = R.value;
      E[I] = SN.value * D[I + 1];
      D[I + 1] = CS.value * D[I + 1];
      if (ICOMPQ == 1) {
        Q[I + 2 * N] = CS.value;
        Q[I + 3 * N] = SN.value;
      } else if (ICOMPQ == 2) {
        WORK[I] = CS.value;
        WORK[NM1 + I] = -SN.value;
      }
    }
  }

  // If ICOMPQ = 0, use DLASDQ to compute the singular values.

  if (ICOMPQ == 0) {
    // Ignore WSTART, instead using WORK[ 1 ], since the two vectors
    // for CS.value and -SN.value above are added only if ICOMPQ == 2,
    // and adding them exceeds documented WORK size of 4*n.
    dlasdq('U', 0, N, 0, 0, 0, D, E, VT, LDVT, U, LDU, U, LDU, WORK(1), INFO);
  } else

  // If N is smaller than the minimum divide size SMLSIZ, then solve
  // the problem with another solver.

  if (N <= SMLSIZ) {
    if (ICOMPQ == 2) {
      dlaset('A', N, N, ZERO, ONE, U, LDU);
      dlaset('A', N, N, ZERO, ONE, VT, LDVT);
      dlasdq('U', 0, N, N, N, 0, D, E, VT, LDVT, U, LDU, U, LDU, WORK(WSTART),
          INFO);
    } else if (ICOMPQ == 1) {
      IU = 1;
      IVT = IU + N;
      dlaset('A', N, N, ZERO, ONE, Q(IU + (QSTART - 1) * N).asMatrix(N), N);
      dlaset('A', N, N, ZERO, ONE, Q(IVT + (QSTART - 1) * N).asMatrix(N), N);
      dlasdq(
          'U',
          0,
          N,
          N,
          N,
          0,
          D,
          E,
          Q(IVT + (QSTART - 1) * N).asMatrix(N),
          N,
          Q(IU + (QSTART - 1) * N).asMatrix(N),
          N,
          Q(IU + (QSTART - 1) * N).asMatrix(N),
          N,
          WORK(WSTART),
          INFO);
    }
  } else {
    if (ICOMPQ == 2) {
      dlaset('A', N, N, ZERO, ONE, U, LDU);
      dlaset('A', N, N, ZERO, ONE, VT, LDVT);
    }

    // Scale.

    ORGNRM = dlanst('M', N, D, E);
    if (ORGNRM == ZERO) return;
    dlascl('G', 0, 0, ORGNRM, ONE, N, 1, D.asMatrix(N), N, IERR);
    dlascl('G', 0, 0, ORGNRM, ONE, NM1, 1, E.asMatrix(NM1), NM1, IERR);

    EPS = (0.9) * dlamch('Epsilon');

    MLVL = (log(N / (SMLSIZ + 1)) ~/ log(TWO)) + 1;
    SMLSZP = SMLSIZ + 1;

    if (ICOMPQ == 1) {
      IU = 1;
      IVT = 1 + SMLSIZ;
      DIFL = IVT + SMLSZP;
      DIFR = DIFL + MLVL;
      Z = DIFR + MLVL * 2;
      IC = Z + MLVL;
      IS = IC + 1;
      POLES = IS + 1;
      GIVNUM = POLES + 2 * MLVL;

      K = 1;
      GIVPTR = 2;
      PERM = 3;
      GIVCOL = PERM + MLVL;
    }

    for (I = 1; I <= N; I++) {
      if ((D[I]).abs() < EPS) {
        D[I] = sign(EPS, D[I]).toDouble();
      }
    }

    START = 1;
    SQRE = 0;

    for (I = 1; I <= NM1; I++) {
      if (((E[I]).abs() < EPS) || (I == NM1)) {
        // Subproblem found. First determine its size and then
        // apply divide and conquer on it.

        if (I < NM1) {
          // A subproblem with E[I] small for I < NM1.

          NSIZE = I - START + 1;
        } else if ((E[I]).abs() >= EPS) {
          // A subproblem with E[NM1] not too small but I = NM1.

          NSIZE = N - START + 1;
        } else {
          // A subproblem with E[NM1] small. This implies an
          // 1-by-1 subproblem at D[N]. Solve this 1-by-1 problem
          // first.

          NSIZE = I - START + 1;
          if (ICOMPQ == 2) {
            U[N][N] = sign(ONE, D[N]).toDouble();
            VT[N][N] = ONE;
          } else if (ICOMPQ == 1) {
            Q[N + (QSTART - 1) * N] = sign(ONE, D[N]).toDouble();
            Q[N + (SMLSIZ + QSTART - 1) * N] = ONE;
          }
          D[N] = (D[N]).abs();
        }
        if (ICOMPQ == 2) {
          dlasd0(NSIZE, SQRE, D(START), E(START), U(START, START), LDU,
              VT(START, START), LDVT, SMLSIZ, IWORK, WORK(WSTART), INFO);
        } else {
          dlasda(
              ICOMPQ,
              SMLSIZ,
              NSIZE,
              SQRE,
              D(START),
              E(START),
              Q(START + (IU + QSTART - 2) * N).asMatrix(N),
              N,
              Q(START + (IVT + QSTART - 2) * N).asMatrix(N),
              IQ(START + K * N),
              Q(START + (DIFL + QSTART - 2) * N).asMatrix(N),
              Q(START + (DIFR + QSTART - 2) * N).asMatrix(N),
              Q(START + (Z + QSTART - 2) * N).asMatrix(N),
              Q(START + (POLES + QSTART - 2) * N).asMatrix(N),
              IQ(START + GIVPTR * N),
              IQ(START + GIVCOL * N).asMatrix(N),
              N,
              IQ(START + PERM * N).asMatrix(N),
              Q(START + (GIVNUM + QSTART - 2) * N).asMatrix(N),
              Q(START + (IC + QSTART - 2) * N),
              Q(START + (IS + QSTART - 2) * N),
              WORK(WSTART),
              IWORK,
              INFO);
        }
        if (INFO.value != 0) {
          return;
        }
        START = I + 1;
      }
    }

    // Unscale
    dlascl('G', 0, 0, ONE, ORGNRM, N, 1, D.asMatrix(N), N, IERR);
  }

  // Use Selection Sort to minimize swaps of singular vectors

  for (II = 2; II <= N; II++) {
    I = II - 1;
    KK = I;
    P = D[I];
    for (J = II; J <= N; J++) {
      if (D[J] > P) {
        KK = J;
        P = D[J];
      }
    }
    if (KK != I) {
      D[KK] = D[I];
      D[I] = P;
      if (ICOMPQ == 1) {
        IQ[I] = KK;
      } else if (ICOMPQ == 2) {
        dswap(N, U(1, I).asArray(), 1, U(1, KK).asArray(), 1);
        dswap(N, VT(I, 1).asArray(), LDVT, VT(KK, 1).asArray(), LDVT);
      }
    } else if (ICOMPQ == 1) {
      IQ[I] = I;
    }
  }

  // If ICOMPQ = 1, use IQ[N,1] as the indicator for UPLO

  if (ICOMPQ == 1) {
    if (IUPLO == 1) {
      IQ[N] = 1;
    } else {
      IQ[N] = 0;
    }
  }

  // If B is lower bidiagonal, update U by those Givens rotations
  // which rotated B to be upper bidiagonal

  if ((IUPLO == 2) && (ICOMPQ == 2)) {
    dlasr('L', 'V', 'B', N, N, WORK(1), WORK(N), U, LDU);
  }
}
