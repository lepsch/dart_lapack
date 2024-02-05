import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgehrd.dart';
import 'package:lapack/src/dhseqr.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dtrevc.dart';
import 'package:lapack/src/dtrsna.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dget37(
  final Array<double> RMAX,
  final Array<int> LMAX,
  final Array<int> NINFO,
  final Box<int> KNT,
  final int NIN,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  const EPSIN = 5.9605e-8;
  const LDT = 20, LWORK = 2 * LDT * (10 + LDT);
  int I, ICMP, IFND, ISCL, J, KMIN, M = 0, N = 0;
  double BIGNUM, EPS, SMLNUM, TNRM, TOL, TOLIN, V, VIMIN, VMAX, VMUL, VRMIN;
  final SELECT = Array<bool>(LDT);
  final IWORK = Array<int>(2 * LDT), LCMP = Array<int>(3);
  final LE = Matrix<double>(LDT, LDT),
      RE = Matrix<double>(LDT, LDT),
      T = Matrix<double>(LDT, LDT),
      TMP = Matrix<double>(LDT, LDT);

  final DUM = Array<double>(1),
      S = Array<double>(LDT),
      SEP = Array<double>(LDT),
      SEPIN = Array<double>(LDT),
      SEPTMP = Array<double>(LDT),
      SIN = Array<double>(LDT),
      STMP = Array<double>(LDT),
      VAL = Array<double>(3),
      WI = Array<double>(LDT),
      WIIN = Array<double>(LDT),
      WITMP = Array<double>(LDT),
      WORK = Array<double>(LWORK),
      WR = Array<double>(LDT),
      WRIN = Array<double>(LDT),
      WRTMP = Array<double>(LDT);
  final INFO = Box(0);

  EPS = dlamch('P');
  SMLNUM = dlamch('S') / EPS;
  BIGNUM = ONE / SMLNUM;

  // EPSIN = 2**(-24) = precision to which input data computed

  EPS = max(EPS, EPSIN);
  RMAX[1] = ZERO;
  RMAX[2] = ZERO;
  RMAX[3] = ZERO;
  LMAX[1] = 0;
  LMAX[2] = 0;
  LMAX[3] = 0;
  KNT.value = 0;
  NINFO[1] = 0;
  NINFO[2] = 0;
  NINFO[3] = 0;

  VAL[1] = sqrt(SMLNUM);
  VAL[2] = ONE;
  VAL[3] = sqrt(BIGNUM);

  // Read input data until N=0.  Assume input eigenvalues are sorted
  // lexicographically (increasing by real part, then decreasing by
  // imaginary part)

  while (true) {
    READ( NIN, FMT = * )N;
    if (N == 0) return;
    for (I = 1; I <= N; I++) {
       READ( NIN, FMT = * )( TMP[ I][ J ], J = 1, N );
    }
    for (I = 1; I <= N; I++) {
       READ( NIN, FMT = * )WRIN[ I ], WIIN[ I ], SIN[ I ], SEPIN[ I ];
    }
    TNRM = dlange('M', N, N, TMP, LDT, WORK);

    // Begin test
    for (ISCL = 1; ISCL <= 3; ISCL++) {
      // Scale input matrix

      KNT.value = KNT.value + 1;
      dlacpy('F', N, N, TMP, LDT, T, LDT);
      VMUL = VAL[ISCL];
      for (I = 1; I <= N; I++) {
        dscal(N, VMUL, T(1, I).asArray(), 1);
      }
      if (TNRM == ZERO) VMUL = ONE;

      // Compute eigenvalues and eigenvectors

      dgehrd(N, 1, N, T, LDT, WORK(1), WORK(N + 1), LWORK - N, INFO);
      if (INFO.value != 0) {
        LMAX[1] = KNT.value;
        NINFO[1] = NINFO[1] + 1;
        continue;
      }
      for (J = 1; J <= N - 2; J++) {
        for (I = J + 2; I <= N; I++) {
          T[I][J] = ZERO;
        }
      }

      // Compute Schur form

      dhseqr(
        'S',
        'N',
        N,
        1,
        N,
        T,
        LDT,
        WR,
        WI,
        DUM,
        1,
        WORK,
        LWORK,
        INFO,
      );
      if (INFO.value != 0) {
        LMAX[2] = KNT.value;
        NINFO[2] = NINFO[2] + 1;
        continue;
      }

      // Compute eigenvectors

      dtrevc(
        'Both',
        'All',
        SELECT,
        N,
        T,
        LDT,
        LE,
        LDT,
        RE,
        LDT,
        N,
        M,
        WORK,
        INFO,
      );

      // Compute condition numbers

      dtrsna(
        'Both',
        'All',
        SELECT,
        N,
        T,
        LDT,
        LE,
        LDT,
        RE,
        LDT,
        S,
        SEP,
        N,
        M,
        WORK,
        N,
        IWORK,
        INFO,
      );
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3] = NINFO[3] + 1;
        continue;
      }

      // Sort eigenvalues and condition numbers lexicographically
      // to compare with inputs

      dcopy(N, WR, 1, WRTMP, 1);
      dcopy(N, WI, 1, WITMP, 1);
      dcopy(N, S, 1, STMP, 1);
      dcopy(N, SEP, 1, SEPTMP, 1);
      dscal(N, ONE / VMUL, SEPTMP, 1);
      for (I = 1; I <= N - 1; I++) {
        KMIN = I;
        VRMIN = WRTMP[I];
        VIMIN = WITMP[I];
        for (J = I + 1; J <= N; J++) {
          if (WRTMP[J] < VRMIN) {
            KMIN = J;
            VRMIN = WRTMP[J];
            VIMIN = WITMP[J];
          }
        }
        WRTMP[KMIN] = WRTMP[I];
        WITMP[KMIN] = WITMP[I];
        WRTMP[I] = VRMIN;
        WITMP[I] = VIMIN;
        VRMIN = STMP[KMIN];
        STMP[KMIN] = STMP[I];
        STMP[I] = VRMIN;
        VRMIN = SEPTMP[KMIN];
        SEPTMP[KMIN] = SEPTMP[I];
        SEPTMP[I] = VRMIN;
      }

      // Compare condition numbers for eigenvalues
      // taking their condition numbers into account

      V = max(TWO * N.toDouble() * EPS * TNRM, SMLNUM);
      if (TNRM == ZERO) V = ONE;
      for (I = 1; I <= N; I++) {
        if (V > SEPTMP[I]) {
          TOL = ONE;
        } else {
          TOL = V / SEPTMP[I];
        }
        if (V > SEPIN[I]) {
          TOLIN = ONE;
        } else {
          TOLIN = V / SEPIN[I];
        }
        TOL = max(TOL, SMLNUM / EPS);
        TOLIN = max(TOLIN, SMLNUM / EPS);
        if (EPS * (SIN[I] - TOLIN) > STMP[I] + TOL) {
          VMAX = ONE / EPS;
        } else if (SIN[I] - TOLIN > STMP[I] + TOL) {
          VMAX = (SIN[I] - TOLIN) / (STMP[I] + TOL);
        } else if (SIN[I] + TOLIN < EPS * (STMP[I] - TOL)) {
          VMAX = ONE / EPS;
        } else if (SIN[I] + TOLIN < STMP[I] - TOL) {
          VMAX = (STMP[I] - TOL) / (SIN[I] + TOLIN);
        } else {
          VMAX = ONE;
        }
        if (VMAX > RMAX[2]) {
          RMAX[2] = VMAX;
          if (NINFO[2] == 0) LMAX[2] = KNT.value;
        }
      }

      // Compare condition numbers for eigenvectors
      // taking their condition numbers into account

      for (I = 1; I <= N; I++) {
        if (V > SEPTMP[I] * STMP[I]) {
          TOL = SEPTMP[I];
        } else {
          TOL = V / STMP[I];
        }
        if (V > SEPIN[I] * SIN[I]) {
          TOLIN = SEPIN[I];
        } else {
          TOLIN = V / SIN[I];
        }
        TOL = max(TOL, SMLNUM / EPS);
        TOLIN = max(TOLIN, SMLNUM / EPS);
        if (EPS * (SEPIN[I] - TOLIN) > SEPTMP[I] + TOL) {
          VMAX = ONE / EPS;
        } else if (SEPIN[I] - TOLIN > SEPTMP[I] + TOL) {
          VMAX = (SEPIN[I] - TOLIN) / (SEPTMP[I] + TOL);
        } else if (SEPIN[I] + TOLIN < EPS * (SEPTMP[I] - TOL)) {
          VMAX = ONE / EPS;
        } else if (SEPIN[I] + TOLIN < SEPTMP[I] - TOL) {
          VMAX = (SEPTMP[I] - TOL) / (SEPIN[I] + TOLIN);
        } else {
          VMAX = ONE;
        }
        if (VMAX > RMAX[2]) {
          RMAX[2] = VMAX;
          if (NINFO[2] == 0) LMAX[2] = KNT.value;
        }
      }

      // Compare condition numbers for eigenvalues
      // without taking their condition numbers into account

      for (I = 1; I <= N; I++) {
        if (SIN[I] <= (2 * N).toDouble() * EPS &&
            STMP[I] <= (2 * N).toDouble() * EPS) {
          VMAX = ONE;
        } else if (EPS * SIN[I] > STMP[I]) {
          VMAX = ONE / EPS;
        } else if (SIN[I] > STMP[I]) {
          VMAX = SIN[I] / STMP[I];
        } else if (SIN[I] < EPS * STMP[I]) {
          VMAX = ONE / EPS;
        } else if (SIN[I] < STMP[I]) {
          VMAX = STMP[I] / SIN[I];
        } else {
          VMAX = ONE;
        }
        if (VMAX > RMAX[3]) {
          RMAX[3] = VMAX;
          if (NINFO[3] == 0) LMAX[3] = KNT.value;
        }
      }

      // Compare condition numbers for eigenvectors
      // without taking their condition numbers into account

      for (I = 1; I <= N; I++) {
        if (SEPIN[I] <= V && SEPTMP[I] <= V) {
          VMAX = ONE;
        } else if (EPS * SEPIN[I] > SEPTMP[I]) {
          VMAX = ONE / EPS;
        } else if (SEPIN[I] > SEPTMP[I]) {
          VMAX = SEPIN[I] / SEPTMP[I];
        } else if (SEPIN[I] < EPS * SEPTMP[I]) {
          VMAX = ONE / EPS;
        } else if (SEPIN[I] < SEPTMP[I]) {
          VMAX = SEPTMP[I] / SEPIN[I];
        } else {
          VMAX = ONE;
        }
        if (VMAX > RMAX[3]) {
          RMAX[3] = VMAX;
          if (NINFO[3] == 0) LMAX[3] = KNT.value;
        }
      }

      // Compute eigenvalue condition numbers only and compare

      VMAX = ZERO;
      DUM[1] = -ONE;
      dcopy(N, DUM, 0, STMP, 1);
      dcopy(N, DUM, 0, SEPTMP, 1);
      dtrsna(
        'Eigcond',
        'All',
        SELECT,
        N,
        T,
        LDT,
        LE,
        LDT,
        RE,
        LDT,
        STMP,
        SEPTMP,
        N,
        M,
        WORK,
        N,
        IWORK,
        INFO.value,
      );
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3] = NINFO[3] + 1;
        continue;
      }
      for (I = 1; I <= N; I++) {
        if (STMP[I] != S[I]) VMAX = ONE / EPS;
        if (SEPTMP[I] != DUM[1]) VMAX = ONE / EPS;
      }

      // Compute eigenvector condition numbers only and compare

      dcopy(N, DUM, 0, STMP, 1);
      dcopy(N, DUM, 0, SEPTMP, 1);
      dtrsna(
        'Veccond',
        'All',
        SELECT,
        N,
        T,
        LDT,
        LE,
        LDT,
        RE,
        LDT,
        STMP,
        SEPTMP,
        N,
        M,
        WORK,
        N,
        IWORK,
        INFO.value,
      );
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3] = NINFO[3] + 1;
        continue;
      }
      for (I = 1; I <= N; I++) {
        if (STMP[I] != DUM[1]) VMAX = ONE / EPS;
        if (SEPTMP[I] != SEP[I]) VMAX = ONE / EPS;
      }

      // Compute all condition numbers using SELECT and compare

      for (I = 1; I <= N; I++) {
        SELECT[I] = true;
      }
      dcopy(N, DUM, 0, STMP, 1);
      dcopy(N, DUM, 0, SEPTMP, 1);
      dtrsna(
        'Bothcond',
        'Some',
        SELECT,
        N,
        T,
        LDT,
        LE,
        LDT,
        RE,
        LDT,
        STMP,
        SEPTMP,
        N,
        M,
        WORK,
        N,
        IWORK,
        INFO.value,
      );
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3] = NINFO[3] + 1;
        continue;
      }
      for (I = 1; I <= N; I++) {
        if (SEPTMP[I] != SEP[I]) VMAX = ONE / EPS;
        if (STMP[I] != S[I]) VMAX = ONE / EPS;
      }

      // Compute eigenvalue condition numbers using SELECT and compare

      dcopy(N, DUM, 0, STMP, 1);
      dcopy(N, DUM, 0, SEPTMP, 1);
      dtrsna(
        'Eigcond',
        'Some',
        SELECT,
        N,
        T,
        LDT,
        LE,
        LDT,
        RE,
        LDT,
        STMP,
        SEPTMP,
        N,
        M,
        WORK,
        N,
        IWORK,
        INFO.value,
      );
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3] = NINFO[3] + 1;
        continue;
      }
      for (I = 1; I <= N; I++) {
        if (STMP[I] != S[I]) VMAX = ONE / EPS;
        if (SEPTMP[I] != DUM[1]) VMAX = ONE / EPS;
      }

      // Compute eigenvector condition numbers using SELECT and compare

      dcopy(N, DUM, 0, STMP, 1);
      dcopy(N, DUM, 0, SEPTMP, 1);
      dtrsna(
        'Veccond',
        'Some',
        SELECT,
        N,
        T,
        LDT,
        LE,
        LDT,
        RE,
        LDT,
        STMP,
        SEPTMP,
        N,
        M,
        WORK,
        N,
        IWORK,
        INFO.value,
      );
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3] = NINFO[3] + 1;
        continue;
      }
      for (I = 1; I <= N; I++) {
        if (STMP[I] != DUM[1]) VMAX = ONE / EPS;
        if (SEPTMP[I] != SEP[I]) VMAX = ONE / EPS;
      }
      if (VMAX > RMAX[1]) {
        RMAX[1] = VMAX;
        if (NINFO[1] == 0) LMAX[1] = KNT.value;
      }

      // Select first real and first complex eigenvalue

      if (WI[1] == ZERO) {
        LCMP[1] = 1;
        IFND = 0;
        for (I = 2; I <= N; I++) {
          if (IFND == 1 || WI[I] == ZERO) {
            SELECT[I] = false;
          } else {
            IFND = 1;
            LCMP[2] = I;
            LCMP[3] = I + 1;
            dcopy(N, RE(1, I).asArray(), 1, RE(1, 2).asArray(), 1);
            dcopy(N, RE(1, I + 1).asArray(), 1, RE(1, 3).asArray(), 1);
            dcopy(N, LE(1, I).asArray(), 1, LE(1, 2).asArray(), 1);
            dcopy(N, LE(1, I + 1).asArray(), 1, LE(1, 3).asArray(), 1);
          }
        }
        if (IFND == 0) {
          ICMP = 1;
        } else {
          ICMP = 3;
        }
      } else {
        LCMP[1] = 1;
        LCMP[2] = 2;
        IFND = 0;
        for (I = 3; I <= N; I++) {
          if (IFND == 1 || WI[I] != ZERO) {
            SELECT[I] = false;
          } else {
            LCMP[3] = I;
            IFND = 1;
            dcopy(N, RE(1, I).asArray(), 1, RE(1, 3).asArray(), 1);
            dcopy(N, LE(1, I).asArray(), 1, LE(1, 3).asArray(), 1);
          }
        }
        if (IFND == 0) {
          ICMP = 2;
        } else {
          ICMP = 3;
        }
      }

      // Compute all selected condition numbers

      dcopy(ICMP, DUM, 0, STMP, 1);
      dcopy(ICMP, DUM, 0, SEPTMP, 1);
      dtrsna(
        'Bothcond',
        'Some',
        SELECT,
        N,
        T,
        LDT,
        LE,
        LDT,
        RE,
        LDT,
        STMP,
        SEPTMP,
        N,
        M,
        WORK,
        N,
        IWORK,
        INFO.value,
      );
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3] = NINFO[3] + 1;
        continue;
      }
      for (I = 1; I <= ICMP; I++) {
        J = LCMP[I];
        if (SEPTMP[I] != SEP[J]) VMAX = ONE / EPS;
        if (STMP[I] != S[J]) VMAX = ONE / EPS;
      }

      // Compute selected eigenvalue condition numbers

      dcopy(ICMP, DUM, 0, STMP, 1);
      dcopy(ICMP, DUM, 0, SEPTMP, 1);
      dtrsna(
        'Eigcond',
        'Some',
        SELECT,
        N,
        T,
        LDT,
        LE,
        LDT,
        RE,
        LDT,
        STMP,
        SEPTMP,
        N,
        M,
        WORK,
        N,
        IWORK,
        INFO.value,
      );
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3] = NINFO[3] + 1;
        continue;
      }
      for (I = 1; I <= ICMP; I++) {
        J = LCMP[I];
        if (STMP[I] != S[J]) VMAX = ONE / EPS;
        if (SEPTMP[I] != DUM[1]) VMAX = ONE / EPS;
      }

      // Compute selected eigenvector condition numbers

      dcopy(ICMP, DUM, 0, STMP, 1);
      dcopy(ICMP, DUM, 0, SEPTMP, 1);
      dtrsna(
        'Veccond',
        'Some',
        SELECT,
        N,
        T,
        LDT,
        LE,
        LDT,
        RE,
        LDT,
        STMP,
        SEPTMP,
        N,
        M,
        WORK,
        N,
        IWORK,
        INFO.value,
      );
      if (INFO.value != 0) {
        LMAX[3] = KNT.value;
        NINFO[3] = NINFO[3] + 1;
        continue;
      }
      for (I = 1; I <= ICMP; I++) {
        J = LCMP[I];
        if (STMP[I] != DUM[1]) VMAX = ONE / EPS;
        if (SEPTMP[I] != SEP[J]) VMAX = ONE / EPS;
      }
      if (VMAX > RMAX[1]) {
        RMAX[1] = VMAX;
        if (NINFO[1] == 0) LMAX[1] = KNT.value;
      }
    }
  }
}
