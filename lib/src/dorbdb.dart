// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/daxpy.dart';
import 'package:dart_lapack/src/blas/dnrm2.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlarf.dart';
import 'package:dart_lapack/src/dlarfgp.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dorbdb(
  final String TRANS,
  final String SIGNS,
  final int M,
  final int P,
  final int Q,
  final Matrix<double> X11_,
  final int LDX11,
  final Matrix<double> X12_,
  final int LDX12,
  final Matrix<double> X21_,
  final int LDX21,
  final Matrix<double> X22_,
  final int LDX22,
  final Array<double> THETA_,
  final Array<double> PHI_,
  final Array<double> TAUP1_,
  final Array<double> TAUP2_,
  final Array<double> TAUQ1_,
  final Array<double> TAUQ2_,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final X11 = X11_.having(ld: LDX11);
  final X12 = X12_.having(ld: LDX12);
  final X21 = X21_.having(ld: LDX21);
  final X22 = X22_.having(ld: LDX22);
  final THETA = THETA_.having();
  final PHI = PHI_.having();
  final TAUP1 = TAUP1_.having();
  final TAUP2 = TAUP2_.having();
  final TAUQ1 = TAUQ1_.having();
  final TAUQ2 = TAUQ2_.having();
  final WORK = WORK_.having();
  const REALONE = 1.0;
  const ONE = 1.0;
  bool COLMAJOR, LQUERY;
  int I, LWORKMIN, LWORKOPT;
  double Z1, Z2, Z3, Z4;

  // Test input arguments

  INFO.value = 0;
  COLMAJOR = !lsame(TRANS, 'T');
  if (!lsame(SIGNS, 'O')) {
    Z1 = REALONE;
    Z2 = REALONE;
    Z3 = REALONE;
    Z4 = REALONE;
  } else {
    Z1 = REALONE;
    Z2 = -REALONE;
    Z3 = REALONE;
    Z4 = -REALONE;
  }
  LQUERY = LWORK == -1;

  if (M < 0) {
    INFO.value = -3;
  } else if (P < 0 || P > M) {
    INFO.value = -4;
  } else if (Q < 0 || Q > P || Q > M - P || Q > M - Q) {
    INFO.value = -5;
  } else if (COLMAJOR && LDX11 < max(1, P)) {
    INFO.value = -7;
  } else if (!COLMAJOR && LDX11 < max(1, Q)) {
    INFO.value = -7;
  } else if (COLMAJOR && LDX12 < max(1, P)) {
    INFO.value = -9;
  } else if (!COLMAJOR && LDX12 < max(1, M - Q)) {
    INFO.value = -9;
  } else if (COLMAJOR && LDX21 < max(1, M - P)) {
    INFO.value = -11;
  } else if (!COLMAJOR && LDX21 < max(1, Q)) {
    INFO.value = -11;
  } else if (COLMAJOR && LDX22 < max(1, M - P)) {
    INFO.value = -13;
  } else if (!COLMAJOR && LDX22 < max(1, M - Q)) {
    INFO.value = -13;
  }

  // Compute workspace

  if (INFO.value == 0) {
    LWORKOPT = M - Q;
    LWORKMIN = M - Q;
    WORK[1] = LWORKOPT.toDouble();
    if (LWORK < LWORKMIN && !LQUERY) {
      INFO.value = -21;
    }
  }
  if (INFO.value != 0) {
    xerbla('xORBDB', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Handle column-major and row-major separately

  if (COLMAJOR) {
    // Reduce columns 1, ..., Q of X11, X12, X21, and X22

    for (I = 1; I <= Q; I++) {
      if (I == 1) {
        dscal(P - I + 1, Z1, X11(I, I).asArray(), 1);
      } else {
        dscal(P - I + 1, Z1 * cos(PHI[I - 1]), X11(I, I).asArray(), 1);
        daxpy(P - I + 1, -Z1 * Z3 * Z4 * sin(PHI[I - 1]),
            X12(I, I - 1).asArray(), 1, X11(I, I).asArray(), 1);
      }
      if (I == 1) {
        dscal(M - P - I + 1, Z2, X21(I, I).asArray(), 1);
      } else {
        dscal(M - P - I + 1, Z2 * cos(PHI[I - 1]), X21(I, I).asArray(), 1);
        daxpy(M - P - I + 1, -Z2 * Z3 * Z4 * sin(PHI[I - 1]),
            X22(I, I - 1).asArray(), 1, X21(I, I).asArray(), 1);
      }

      THETA[I] = atan2(dnrm2(M - P - I + 1, X21(I, I).asArray(), 1),
          dnrm2(P - I + 1, X11(I, I).asArray(), 1));

      if (P > I) {
        dlarfgp(
            P - I + 1, X11.box(I, I), X11(I + 1, I).asArray(), 1, TAUP1.box(I));
      } else if (P == I) {
        dlarfgp(P - I + 1, X11.box(I, I), X11(I, I).asArray(), 1, TAUP1.box(I));
      }
      X11[I][I] = ONE;
      if (M - P > I) {
        dlarfgp(M - P - I + 1, X21.box(I, I), X21(I + 1, I).asArray(), 1,
            TAUP2.box(I));
      } else if (M - P == I) {
        dlarfgp(
            M - P - I + 1, X21.box(I, I), X21(I, I).asArray(), 1, TAUP2.box(I));
      }
      X21[I][I] = ONE;

      if (Q > I) {
        dlarf('L', P - I + 1, Q - I, X11(I, I).asArray(), 1, TAUP1[I],
            X11(I, I + 1), LDX11, WORK);
      }
      if (M - Q + 1 > I) {
        dlarf('L', P - I + 1, M - Q - I + 1, X11(I, I).asArray(), 1, TAUP1[I],
            X12(I, I), LDX12, WORK);
      }
      if (Q > I) {
        dlarf('L', M - P - I + 1, Q - I, X21(I, I).asArray(), 1, TAUP2[I],
            X21(I, I + 1), LDX21, WORK);
      }
      if (M - Q + 1 > I) {
        dlarf('L', M - P - I + 1, M - Q - I + 1, X21(I, I).asArray(), 1,
            TAUP2[I], X22(I, I), LDX22, WORK);
      }

      if (I < Q) {
        dscal(Q - I, -Z1 * Z3 * sin(THETA[I]), X11(I, I + 1).asArray(), LDX11);
        daxpy(Q - I, Z2 * Z3 * cos(THETA[I]), X21(I, I + 1).asArray(), LDX21,
            X11(I, I + 1).asArray(), LDX11);
      }
      dscal(
          M - Q - I + 1, -Z1 * Z4 * sin(THETA[I]), X12(I, I).asArray(), LDX12);
      daxpy(M - Q - I + 1, Z2 * Z4 * cos(THETA[I]), X22(I, I).asArray(), LDX22,
          X12(I, I).asArray(), LDX12);

      if (I < Q) {
        PHI[I] = atan2(dnrm2(Q - I, X11(I, I + 1).asArray(), LDX11),
            dnrm2(M - Q - I + 1, X12(I, I).asArray(), LDX12));
      }

      if (I < Q) {
        if (Q - I == 1) {
          dlarfgp(Q - I, X11.box(I, I + 1), X11(I, I + 1).asArray(), LDX11,
              TAUQ1.box(I));
        } else {
          dlarfgp(Q - I, X11.box(I, I + 1), X11(I, I + 2).asArray(), LDX11,
              TAUQ1.box(I));
        }
        X11[I][I + 1] = ONE;
      }
      if (Q + I - 1 < M) {
        if (M - Q == I) {
          dlarfgp(M - Q - I + 1, X12.box(I, I), X12(I, I).asArray(), LDX12,
              TAUQ2.box(I));
        } else {
          dlarfgp(M - Q - I + 1, X12.box(I, I), X12(I, I + 1).asArray(), LDX12,
              TAUQ2.box(I));
        }
      }
      X12[I][I] = ONE;

      if (I < Q) {
        dlarf('R', P - I, Q - I, X11(I, I + 1).asArray(), LDX11, TAUQ1[I],
            X11(I + 1, I + 1), LDX11, WORK);
        dlarf('R', M - P - I, Q - I, X11(I, I + 1).asArray(), LDX11, TAUQ1[I],
            X21(I + 1, I + 1), LDX21, WORK);
      }
      if (P > I) {
        dlarf('R', P - I, M - Q - I + 1, X12(I, I).asArray(), LDX12, TAUQ2[I],
            X12(I + 1, I), LDX12, WORK);
      }
      if (M - P > I) {
        dlarf('R', M - P - I, M - Q - I + 1, X12(I, I).asArray(), LDX12,
            TAUQ2[I], X22(I + 1, I), LDX22, WORK);
      }
    }

    // Reduce columns Q + 1, ..., P of X12, X22

    for (I = Q + 1; I <= P; I++) {
      dscal(M - Q - I + 1, -Z1 * Z4, X12(I, I).asArray(), LDX12);
      if (I >= M - Q) {
        dlarfgp(M - Q - I + 1, X12.box(I, I), X12(I, I).asArray(), LDX12,
            TAUQ2.box(I));
      } else {
        dlarfgp(M - Q - I + 1, X12.box(I, I), X12(I, I + 1).asArray(), LDX12,
            TAUQ2.box(I));
      }
      X12[I][I] = ONE;

      if (P > I) {
        dlarf('R', P - I, M - Q - I + 1, X12(I, I).asArray(), LDX12, TAUQ2[I],
            X12(I + 1, I), LDX12, WORK);
      }
      if (M - P - Q >= 1) {
        dlarf('R', M - P - Q, M - Q - I + 1, X12(I, I).asArray(), LDX12,
            TAUQ2[I], X22(Q + 1, I), LDX22, WORK);
      }
    }

    // Reduce columns P + 1, ..., M - Q of X12, X22

    for (I = 1; I <= M - P - Q; I++) {
      dscal(M - P - Q - I + 1, Z2 * Z4, X22(Q + I, P + I).asArray(), LDX22);
      if (I == M - P - Q) {
        dlarfgp(M - P - Q - I + 1, X22.box(Q + I, P + I),
            X22(Q + I, P + I).asArray(), LDX22, TAUQ2.box(P + I));
      } else {
        dlarfgp(M - P - Q - I + 1, X22.box(Q + I, P + I),
            X22(Q + I, P + I + 1).asArray(), LDX22, TAUQ2.box(P + I));
      }
      X22[Q + I][P + I] = ONE;
      if (I < M - P - Q) {
        dlarf(
            'R',
            M - P - Q - I,
            M - P - Q - I + 1,
            X22(Q + I, P + I).asArray(),
            LDX22,
            TAUQ2[P + I],
            X22(Q + I + 1, P + I),
            LDX22,
            WORK);
      }
    }
  } else {
    // Reduce columns 1, ..., Q of X11, X12, X21, X22

    for (I = 1; I <= Q; I++) {
      if (I == 1) {
        dscal(P - I + 1, Z1, X11(I, I).asArray(), LDX11);
      } else {
        dscal(P - I + 1, Z1 * cos(PHI[I - 1]), X11(I, I).asArray(), LDX11);
        daxpy(P - I + 1, -Z1 * Z3 * Z4 * sin(PHI[I - 1]),
            X12(I - 1, I).asArray(), LDX12, X11(I, I).asArray(), LDX11);
      }
      if (I == 1) {
        dscal(M - P - I + 1, Z2, X21(I, I).asArray(), LDX21);
      } else {
        dscal(M - P - I + 1, Z2 * cos(PHI[I - 1]), X21(I, I).asArray(), LDX21);
        daxpy(M - P - I + 1, -Z2 * Z3 * Z4 * sin(PHI[I - 1]),
            X22(I - 1, I).asArray(), LDX22, X21(I, I).asArray(), LDX21);
      }

      THETA[I] = atan2(dnrm2(M - P - I + 1, X21(I, I).asArray(), LDX21),
          dnrm2(P - I + 1, X11(I, I).asArray(), LDX11));

      dlarfgp(P - I + 1, X11.box(I, I), X11(I, I + 1).asArray(), LDX11,
          TAUP1.box(I));
      X11[I][I] = ONE;
      if (I == M - P) {
        dlarfgp(M - P - I + 1, X21.box(I, I), X21(I, I).asArray(), LDX21,
            TAUP2.box(I));
      } else {
        dlarfgp(M - P - I + 1, X21.box(I, I), X21(I, I + 1).asArray(), LDX21,
            TAUP2.box(I));
      }
      X21[I][I] = ONE;

      if (Q > I) {
        dlarf('R', Q - I, P - I + 1, X11(I, I).asArray(), LDX11, TAUP1[I],
            X11(I + 1, I), LDX11, WORK);
      }
      if (M - Q + 1 > I) {
        dlarf('R', M - Q - I + 1, P - I + 1, X11(I, I).asArray(), LDX11,
            TAUP1[I], X12(I, I), LDX12, WORK);
      }
      if (Q > I) {
        dlarf('R', Q - I, M - P - I + 1, X21(I, I).asArray(), LDX21, TAUP2[I],
            X21(I + 1, I), LDX21, WORK);
      }
      if (M - Q + 1 > I) {
        dlarf('R', M - Q - I + 1, M - P - I + 1, X21(I, I).asArray(), LDX21,
            TAUP2[I], X22(I, I), LDX22, WORK);
      }

      if (I < Q) {
        dscal(Q - I, -Z1 * Z3 * sin(THETA[I]), X11(I + 1, I).asArray(), 1);
        daxpy(Q - I, Z2 * Z3 * cos(THETA[I]), X21(I + 1, I).asArray(), 1,
            X11(I + 1, I).asArray(), 1);
      }
      dscal(M - Q - I + 1, -Z1 * Z4 * sin(THETA[I]), X12(I, I).asArray(), 1);
      daxpy(M - Q - I + 1, Z2 * Z4 * cos(THETA[I]), X22(I, I).asArray(), 1,
          X12(I, I).asArray(), 1);

      if (I < Q) {
        PHI[I] = atan2(dnrm2(Q - I, X11(I + 1, I).asArray(), 1),
            dnrm2(M - Q - I + 1, X12(I, I).asArray(), 1));
      }

      if (I < Q) {
        if (Q - I == 1) {
          dlarfgp(Q - I, X11.box(I + 1, I), X11(I + 1, I).asArray(), 1,
              TAUQ1.box(I));
        } else {
          dlarfgp(Q - I, X11.box(I + 1, I), X11(I + 2, I).asArray(), 1,
              TAUQ1.box(I));
        }
        X11[I + 1][I] = ONE;
      }
      if (M - Q > I) {
        dlarfgp(M - Q - I + 1, X12.box(I, I), X12(I + 1, I).asArray(), 1,
            TAUQ2.box(I));
      } else {
        dlarfgp(
            M - Q - I + 1, X12.box(I, I), X12(I, I).asArray(), 1, TAUQ2.box(I));
      }
      X12[I][I] = ONE;

      if (I < Q) {
        dlarf('L', Q - I, P - I, X11(I + 1, I).asArray(), 1, TAUQ1[I],
            X11(I + 1, I + 1), LDX11, WORK);
        dlarf('L', Q - I, M - P - I, X11(I + 1, I).asArray(), 1, TAUQ1[I],
            X21(I + 1, I + 1), LDX21, WORK);
      }
      dlarf('L', M - Q - I + 1, P - I, X12(I, I).asArray(), 1, TAUQ2[I],
          X12(I, I + 1), LDX12, WORK);
      if (M - P - I > 0) {
        dlarf('L', M - Q - I + 1, M - P - I, X12(I, I).asArray(), 1, TAUQ2[I],
            X22(I, I + 1), LDX22, WORK);
      }
    }

    // Reduce columns Q + 1, ..., P of X12, X22

    for (I = Q + 1; I <= P; I++) {
      dscal(M - Q - I + 1, -Z1 * Z4, X12(I, I).asArray(), 1);
      dlarfgp(M - Q - I + 1, X12.box(I, I), X12(I + 1, I).asArray(), 1,
          TAUQ2.box(I));
      X12[I][I] = ONE;

      if (P > I) {
        dlarf('L', M - Q - I + 1, P - I, X12(I, I).asArray(), 1, TAUQ2[I],
            X12(I, I + 1), LDX12, WORK);
      }
      if (M - P - Q >= 1) {
        dlarf('L', M - Q - I + 1, M - P - Q, X12(I, I).asArray(), 1, TAUQ2[I],
            X22(I, Q + 1), LDX22, WORK);
      }
    }

    // Reduce columns P + 1, ..., M - Q of X12, X22

    for (I = 1; I <= M - P - Q; I++) {
      dscal(M - P - Q - I + 1, Z2 * Z4, X22(P + I, Q + I).asArray(), 1);
      if (M - P - Q == I) {
        dlarfgp(M - P - Q - I + 1, X22.box(P + I, Q + I),
            X22(P + I, Q + I).asArray(), 1, TAUQ2.box(P + I));
      } else {
        dlarfgp(M - P - Q - I + 1, X22.box(P + I, Q + I),
            X22(P + I + 1, Q + I).asArray(), 1, TAUQ2.box(P + I));
        dlarf(
            'L',
            M - P - Q - I + 1,
            M - P - Q - I,
            X22(P + I, Q + I).asArray(),
            1,
            TAUQ2[P + I],
            X22(P + I, Q + I + 1),
            LDX22,
            WORK);
      }
      X22[P + I][Q + I] = ONE;
    }
  }
}
