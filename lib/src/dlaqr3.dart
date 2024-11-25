// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgehrd.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlahqr.dart';
import 'package:dart_lapack/src/dlanv2.dart';
import 'package:dart_lapack/src/dlaqr4.dart';
import 'package:dart_lapack/src/dlarf.dart';
import 'package:dart_lapack/src/dlarfg.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dormhr.dart';
import 'package:dart_lapack/src/dtrexc.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void dlaqr3(
  final bool WANTT,
  final bool WANTZ,
  final int N,
  final int KTOP,
  final int KBOT,
  final int NW,
  final Matrix<double> H_,
  final int LDH,
  final int ILOZ,
  final int IHIZ,
  final Matrix<double> Z_,
  final int LDZ,
  final Box<int> NS,
  final Box<int> ND,
  final Array<double> SR_,
  final Array<double> SI_,
  final Matrix<double> V_,
  final int LDV,
  final int NH,
  final Matrix<double> T_,
  final int LDT,
  final int NV,
  final Matrix<double> WV_,
  final int LDWV,
  final Array<double> WORK_,
  final int LWORK,
) {
  final H = H_.having(ld: LDH);
  final Z = Z_.having(ld: LDZ);
  final SR = SR_.having();
  final SI = SI_.having();
  final V = V_.having(ld: LDV);
  final T = T_.having(ld: LDT);
  final WV = WV_.having(ld: LDWV);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  double EVI, EVK, FOO, S, SAFMIN, SMLNUM, ULP;
  int I,
      J,
      JW,
      K,
      KCOL,
      KEND,
      KLN,
      KROW,
      KWTOP,
      LTOP,
      LWK1,
      LWK2,
      LWK3,
      LWKOPT,
      NMIN;
  bool BULGE, SORTED;
  final INFO = Box(0), INFQR = Box(0), IFST = Box(0), ILST = Box(0);
  final AA = Box(0.0),
      BB = Box(0.0),
      CC = Box(0.0),
      CS = Box(0.0),
      DD = Box(0.0),
      SN = Box(0.0),
      TAU = Box(0.0),
      BETA = Box(0.0);

  // Estimate optimal workspace.
  JW = min(NW, KBOT - KTOP + 1);
  if (JW <= 2) {
    LWKOPT = 1;
  } else {
    // Workspace query call to DGEHRD
    dgehrd(JW, 1, JW - 1, T, LDT, WORK, WORK, -1, INFO);
    LWK1 = WORK[1].toInt();

    // Workspace query call to DORMHR
    dormhr('R', 'N', JW, JW, 1, JW - 1, T, LDT, WORK, V, LDV, WORK, -1, INFO);
    LWK2 = WORK[1].toInt();

    // Workspace query call to DLAQR4
    dlaqr4(
        true, true, JW, 1, JW, T, LDT, SR, SI, 1, JW, V, LDV, WORK, -1, INFQR);
    LWK3 = WORK[1].toInt();

    // Optimal workspace
    LWKOPT = max(JW + max(LWK1, LWK2), LWK3);
  }

  // Quick return in case of workspace query.
  if (LWORK == -1) {
    WORK[1] = LWKOPT.toDouble();
    return;
  }

  // Nothing to do
  // for an empty active block
  NS.value = 0;
  ND.value = 0;
  WORK[1] = ONE;
  if (KTOP > KBOT) return;
  // nor for an empty deflation window.
  if (NW < 1) return;

  // Machine constants
  SAFMIN = dlamch('SAFE MINIMUM');
  ULP = dlamch('PRECISION');
  SMLNUM = SAFMIN * (N / ULP);

  // Setup deflation window
  JW = min(NW, KBOT - KTOP + 1);
  KWTOP = KBOT - JW + 1;
  if (KWTOP == KTOP) {
    S = ZERO;
  } else {
    S = H[KWTOP][KWTOP - 1];
  }

  if (KBOT == KWTOP) {
    // 1-by-1 deflation window: not much to do
    SR[KWTOP] = H[KWTOP][KWTOP];
    SI[KWTOP] = ZERO;
    NS.value = 1;
    ND.value = 0;
    if (S.abs() <= max(SMLNUM, ULP * H[KWTOP][KWTOP].abs())) {
      NS.value = 0;
      ND.value = 1;
      if (KWTOP > KTOP) H[KWTOP][KWTOP - 1] = ZERO;
    }
    WORK[1] = ONE;
    return;
  }

  // Convert to spike-triangular form.  (In case of a
  // rare QR failure, this routine continues to do
  // aggressive early deflation using that part of
  // the deflation window that converged using INFQR
  // here and there to keep track.)
  dlacpy('U', JW, JW, H(KWTOP, KWTOP), LDH, T, LDT);
  dcopy(JW - 1, H(KWTOP + 1, KWTOP).asArray(), LDH + 1, T(2, 1).asArray(),
      LDT + 1);

  dlaset('A', JW, JW, ZERO, ONE, V, LDV);
  NMIN = ilaenv(12, 'DLAQR3', 'SV', JW, 1, JW, LWORK);
  if (JW > NMIN) {
    dlaqr4(true, true, JW, 1, JW, T, LDT, SR(KWTOP), SI(KWTOP), 1, JW, V, LDV,
        WORK, LWORK, INFQR);
  } else {
    dlahqr(true, true, JW, 1, JW, T, LDT, SR(KWTOP), SI(KWTOP), 1, JW, V, LDV,
        INFQR);
  }

  // DTREXC needs a clean margin near the diagonal
  for (J = 1; J <= JW - 3; J++) {
    T[J + 2][J] = ZERO;
    T[J + 3][J] = ZERO;
  }
  if (JW > 2) T[JW][JW - 2] = ZERO;

  // Deflation detection loop
  NS.value = JW;
  ILST.value = INFQR.value + 1;
  while (ILST.value <= NS.value) {
    if (NS.value == 1) {
      BULGE = false;
    } else {
      BULGE = T[NS.value][NS.value - 1] != ZERO;
    }

    // Small spike tip test for deflation
    if (!BULGE) {
      // Real eigenvalue
      FOO = T[NS.value][NS.value].abs();
      if (FOO == ZERO) FOO = S.abs();
      if ((S * V[1][NS.value]).abs() <= max(SMLNUM, ULP * FOO)) {
        // Deflatable
        NS.value--;
      } else {
        // Undeflatable.   Move it up out of the way.
        // (DTREXC can not fail in this case.)
        IFST.value = NS.value;
        dtrexc('V', JW, T, LDT, V, LDV, IFST, ILST, WORK, INFO);
        ILST.value++;
      }
    } else {
      // Complex conjugate pair
      FOO = T[NS.value][NS.value].abs() +
          sqrt(T[NS.value][NS.value - 1].abs()) *
              sqrt(T[NS.value - 1][NS.value].abs());
      if (FOO == ZERO) FOO = S.abs();
      if (max((S * V[1][NS.value]).abs(), (S * V[1][NS.value - 1]).abs()) <=
          max(SMLNUM, ULP * FOO)) {
        // Deflatable

        NS.value -= 2;
      } else {
        // Undeflatable. Move them up out of the way.
        // Fortunately, DTREXC does the right thing with
        // ILST in case of a rare exchange failure.
        IFST.value = NS.value;
        dtrexc('V', JW, T, LDT, V, LDV, IFST, ILST, WORK, INFO);
        ILST.value += 2;
      }
    }

    // End deflation detection loop
  }

  // Return to Hessenberg form
  if (NS.value == 0) S = ZERO;

  if (NS.value < JW) {
    // sorting diagonal blocks of T improves accuracy for
    // graded matrices.  Bubble sort deals well with
    // exchange failures.
    SORTED = false;
    I = NS.value + 1;
    while (!SORTED) {
      SORTED = true;

      KEND = I - 1;
      I = INFQR.value + 1;
      if (I == NS.value) {
        K = I + 1;
      } else if (T[I + 1][I] == ZERO) {
        K = I + 1;
      } else {
        K = I + 2;
      }
      while (K <= KEND) {
        if (K == I + 1) {
          EVI = T[I][I].abs();
        } else {
          EVI =
              T[I][I].abs() + sqrt(T[I + 1][I].abs()) * sqrt(T[I][I + 1].abs());
        }

        if (K == KEND) {
          EVK = T[K][K].abs();
        } else if (T[K + 1][K] == ZERO) {
          EVK = T[K][K].abs();
        } else {
          EVK =
              T[K][K].abs() + sqrt(T[K + 1][K].abs()) * sqrt(T[K][K + 1].abs());
        }

        if (EVI >= EVK) {
          I = K;
        } else {
          SORTED = false;
          IFST.value = I;
          ILST.value = K;
          dtrexc('V', JW, T, LDT, V, LDV, IFST, ILST, WORK, INFO);
          if (INFO.value == 0) {
            I = ILST.value;
          } else {
            I = K;
          }
        }
        if (I == KEND) {
          K = I + 1;
        } else if (T[I + 1][I] == ZERO) {
          K = I + 1;
        } else {
          K = I + 2;
        }
      }
    }
  }

  // Restore shift/eigenvalue array from T
  I = JW;
  while (I >= INFQR.value + 1) {
    if (I == INFQR.value + 1) {
      SR[KWTOP + I - 1] = T[I][I];
      SI[KWTOP + I - 1] = ZERO;
      I--;
    } else if (T[I][I - 1] == ZERO) {
      SR[KWTOP + I - 1] = T[I][I];
      SI[KWTOP + I - 1] = ZERO;
      I--;
    } else {
      AA.value = T[I - 1][I - 1];
      CC.value = T[I][I - 1];
      BB.value = T[I - 1][I];
      DD.value = T[I][I];
      dlanv2(AA, BB, CC, DD, SR.box(KWTOP + I - 2), SI.box(KWTOP + I - 2),
          SR.box(KWTOP + I - 1), SI.box(KWTOP + I - 1), CS, SN);
      I -= 2;
    }
  }

  if (NS.value < JW || S == ZERO) {
    if (NS.value > 1 && S != ZERO) {
      // Reflect spike back into lower triangle
      dcopy(NS.value, V.asArray(), LDV, WORK, 1);
      BETA.value = WORK[1];
      dlarfg(NS.value, BETA, WORK(2), 1, TAU);
      WORK[1] = ONE;

      dlaset('L', JW - 2, JW - 2, ZERO, ZERO, T(3, 1), LDT);

      dlarf('L', NS.value, JW, WORK, 1, TAU.value, T, LDT, WORK(JW + 1));
      dlarf('R', NS.value, NS.value, WORK, 1, TAU.value, T, LDT, WORK(JW + 1));
      dlarf('R', JW, NS.value, WORK, 1, TAU.value, V, LDV, WORK(JW + 1));

      dgehrd(JW, 1, NS.value, T, LDT, WORK, WORK(JW + 1), LWORK - JW, INFO);
    }

    // Copy updated reduced window into place
    if (KWTOP > 1) H[KWTOP][KWTOP - 1] = S * V[1][1];
    dlacpy('U', JW, JW, T, LDT, H(KWTOP, KWTOP), LDH);
    dcopy(JW - 1, T(2, 1).asArray(), LDT + 1, H(KWTOP + 1, KWTOP).asArray(),
        LDH + 1);

    // Accumulate orthogonal matrix in order update
    // H and Z, if requested.
    if (NS.value > 1 && S != ZERO) {
      dormhr('R', 'N', JW, NS.value, 1, NS.value, T, LDT, WORK, V, LDV,
          WORK(JW + 1), LWORK - JW, INFO);
    }

    // Update vertical slab in H
    if (WANTT) {
      LTOP = 1;
    } else {
      LTOP = KTOP;
    }
    for (KROW = LTOP;
        NV < 0 ? KROW >= KWTOP - 1 : KROW <= KWTOP - 1;
        KROW += NV) {
      KLN = min(NV, KWTOP - KROW);
      dgemm('N', 'N', KLN, JW, JW, ONE, H(KROW, KWTOP), LDH, V, LDV, ZERO, WV,
          LDWV);
      dlacpy('A', KLN, JW, WV, LDWV, H(KROW, KWTOP), LDH);
    }

    // Update horizontal slab in H
    if (WANTT) {
      for (KCOL = KBOT + 1; NH < 0 ? KCOL >= N : KCOL <= N; KCOL += NH) {
        KLN = min(NH, N - KCOL + 1);
        dgemm('C', 'N', JW, KLN, JW, ONE, V, LDV, H(KWTOP, KCOL), LDH, ZERO, T,
            LDT);
        dlacpy('A', JW, KLN, T, LDT, H(KWTOP, KCOL), LDH);
      }
    }

    // Update vertical slab in Z
    if (WANTZ) {
      for (KROW = ILOZ; NV < 0 ? KROW >= IHIZ : KROW <= IHIZ; KROW += NV) {
        KLN = min(NV, IHIZ - KROW + 1);
        dgemm('N', 'N', KLN, JW, JW, ONE, Z(KROW, KWTOP), LDZ, V, LDV, ZERO, WV,
            LDWV);
        dlacpy('A', KLN, JW, WV, LDWV, Z(KROW, KWTOP), LDZ);
      }
    }
  }

  // Return the number of deflations
  ND.value = JW - NS.value;

  // and the number of shifts. (Subtracting
  // INFQR from the spike length takes care
  // of the case of a rare QR failure while
  // calculating eigenvalues of the deflation
  // window.)
  NS.value -= INFQR.value;

  // Return optimal workspace.
  WORK[1] = LWKOPT.toDouble();
}
