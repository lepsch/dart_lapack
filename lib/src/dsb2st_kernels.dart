// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlarfg.dart';
import 'package:dart_lapack/src/dlarfx.dart';
import 'package:dart_lapack/src/dlarfy.dart';
import 'package:dart_lapack/src/matrix.dart';

void dsb2st_kernels(
  final String UPLO,
  final bool WANTZ,
  final int TTYPE,
  final int ST,
  final int ED,
  final int SWEEP,
  final int N,
  final int NB,
  final int IB,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> V_,
  final Array<double> TAU_,
  final int LDVT,
  final Array<double> WORK_,
) {
  final A = A_.having(ld: LDA);
  final V = V_.having();
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool UPPER;
  int I, J1, J2, LM, LN, VPOS, TAUPOS, DPOS, OFDPOS;
  final CTMP = Box(0.0);

  UPPER = lsame(UPLO, 'U');

  if (UPPER) {
    DPOS = 2 * NB + 1;
    OFDPOS = 2 * NB;
  } else {
    DPOS = 1;
    OFDPOS = 2;
  }

  // Upper case

  if (UPPER) {
    if (WANTZ) {
      VPOS = ((SWEEP - 1) % 2) * N + ST;
      TAUPOS = ((SWEEP - 1) % 2) * N + ST;
    } else {
      VPOS = ((SWEEP - 1) % 2) * N + ST;
      TAUPOS = ((SWEEP - 1) % 2) * N + ST;
    }

    if (TTYPE == 1) {
      LM = ED - ST + 1;

      V[VPOS] = ONE;
      for (I = 1; I <= LM - 1; I++) {
        V[VPOS + I] = A[OFDPOS - I][ST + I];
        A[OFDPOS - I][ST + I] = ZERO;
      }
      CTMP.value = A[OFDPOS][ST];
      dlarfg(LM, CTMP, V(VPOS + 1), 1, TAU(TAUPOS));
      A[OFDPOS][ST] = CTMP.value;

      LM = ED - ST + 1;
      dlarfy(UPLO, LM, V(VPOS), 1, TAU[TAUPOS], A(DPOS, ST), LDA - 1, WORK);
    }

    if (TTYPE == 3) {
      LM = ED - ST + 1;
      dlarfy(UPLO, LM, V(VPOS), 1, TAU[TAUPOS], A(DPOS, ST), LDA - 1, WORK);
    }

    if (TTYPE == 2) {
      J1 = ED + 1;
      J2 = min(ED + NB, N);
      LN = ED - ST + 1;
      LM = J2 - J1 + 1;
      if (LM > 0) {
        dlarfx('Left', LN, LM, V(VPOS), TAU[TAUPOS], A(DPOS - NB, J1), LDA - 1,
            WORK);

        if (WANTZ) {
          VPOS = ((SWEEP - 1) % 2) * N + J1;
          TAUPOS = ((SWEEP - 1) % 2) * N + J1;
        } else {
          VPOS = ((SWEEP - 1) % 2) * N + J1;
          TAUPOS = ((SWEEP - 1) % 2) * N + J1;
        }

        V[VPOS] = ONE;
        for (I = 1; I <= LM - 1; I++) {
          V[VPOS + I] = A[DPOS - NB - I][J1 + I];
          A[DPOS - NB - I][J1 + I] = ZERO;
        }
        CTMP.value = A[DPOS - NB][J1];
        dlarfg(LM, CTMP, V(VPOS + 1), 1, TAU(TAUPOS));
        A[DPOS - NB][J1] = CTMP.value;

        dlarfx('Right', LN - 1, LM, V(VPOS), TAU[TAUPOS], A(DPOS - NB + 1, J1),
            LDA - 1, WORK);
      }
    }

    // Lower case
  } else {
    if (WANTZ) {
      VPOS = ((SWEEP - 1) % 2) * N + ST;
      TAUPOS = ((SWEEP - 1) % 2) * N + ST;
    } else {
      VPOS = ((SWEEP - 1) % 2) * N + ST;
      TAUPOS = ((SWEEP - 1) % 2) * N + ST;
    }

    if (TTYPE == 1) {
      LM = ED - ST + 1;

      V[VPOS] = ONE;
      for (I = 1; I <= LM - 1; I++) {
        V[VPOS + I] = A[OFDPOS + I][ST - 1];
        A[OFDPOS + I][ST - 1] = ZERO;
      }
      dlarfg(LM, A(OFDPOS, ST - 1), V(VPOS + 1), 1, TAU(TAUPOS));

      LM = ED - ST + 1;

      dlarfy(UPLO, LM, V(VPOS), 1, TAU[TAUPOS], A(DPOS, ST), LDA - 1, WORK);
    }

    if (TTYPE == 3) {
      LM = ED - ST + 1;

      dlarfy(UPLO, LM, V(VPOS), 1, TAU[TAUPOS], A(DPOS, ST), LDA - 1, WORK);
    }

    if (TTYPE == 2) {
      J1 = ED + 1;
      J2 = min(ED + NB, N);
      LN = ED - ST + 1;
      LM = J2 - J1 + 1;

      if (LM > 0) {
        dlarfx('Right', LM, LN, V(VPOS), TAU[TAUPOS], A(DPOS + NB, ST), LDA - 1,
            WORK);

        if (WANTZ) {
          VPOS = ((SWEEP - 1) % 2) * N + J1;
          TAUPOS = ((SWEEP - 1) % 2) * N + J1;
        } else {
          VPOS = ((SWEEP - 1) % 2) * N + J1;
          TAUPOS = ((SWEEP - 1) % 2) * N + J1;
        }

        V[VPOS] = ONE;
        for (I = 1; I <= LM - 1; I++) {
          V[VPOS + I] = A[DPOS + NB + I][ST];
          A[DPOS + NB + I][ST] = ZERO;
        }
        dlarfg(LM, A(DPOS + NB, ST), V(VPOS + 1), 1, TAU(TAUPOS));

        dlarfx('Left', LM, LN - 1, V(VPOS), TAU[TAUPOS],
            A(DPOS + NB - 1, ST + 1), LDA - 1, WORK);
      }
    }
  }
}
