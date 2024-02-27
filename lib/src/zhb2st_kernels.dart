import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlarfg.dart';
import 'package:lapack/src/zlarfx.dart';
import 'package:lapack/src/zlarfy.dart';

void zhb2st_kernels(
  final String UPLO,
  final bool WANTZ,
  final int TTYPE,
  final int ST,
  final int ED,
  final int SWEEP,
  final int N,
  final int NB,
  final int IB,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> V_,
  final Array<Complex> TAU_,
  final int LDVT,
  final Array<Complex> WORK_,
) {
  final A = A_.dim(LDA);
  final V = V_.dim();
  final TAU = TAU_.dim();
  final WORK = WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  bool UPPER;
  int I, J1, J2, LM, LN, VPOS, TAUPOS, DPOS, OFDPOS;
  final CTMP = Box(Complex.zero);

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

      V[VPOS] = Complex.one;
      for (I = 1; I <= LM - 1; I++) {
        // 10
        V[VPOS + I] = A[OFDPOS - I][ST + I].conjugate();
        A[OFDPOS - I][ST + I] = Complex.zero;
      } // 10
      CTMP.value = A[OFDPOS][ST].conjugate();
      zlarfg(LM, CTMP, V(VPOS + 1), 1, TAU(TAUPOS));
      A[OFDPOS][ST] = CTMP.value;

      LM = ED - ST + 1;
      zlarfy(UPLO, LM, V(VPOS), 1, TAU[TAUPOS].conjugate(), A(DPOS, ST),
          LDA - 1, WORK);
    }

    if (TTYPE == 3) {
      LM = ED - ST + 1;
      zlarfy(UPLO, LM, V(VPOS), 1, TAU[TAUPOS].conjugate(), A(DPOS, ST),
          LDA - 1, WORK);
    }

    if (TTYPE == 2) {
      J1 = ED + 1;
      J2 = min(ED + NB, N);
      LN = ED - ST + 1;
      LM = J2 - J1 + 1;
      if (LM > 0) {
        zlarfx('Left', LN, LM, V(VPOS), TAU[TAUPOS].conjugate(),
            A(DPOS - NB, J1), LDA - 1, WORK);

        if (WANTZ) {
          VPOS = ((SWEEP - 1) % 2) * N + J1;
          TAUPOS = ((SWEEP - 1) % 2) * N + J1;
        } else {
          VPOS = ((SWEEP - 1) % 2) * N + J1;
          TAUPOS = ((SWEEP - 1) % 2) * N + J1;
        }

        V[VPOS] = Complex.one;
        for (I = 1; I <= LM - 1; I++) {
          // 30
          V[VPOS + I] = A[DPOS - NB - I][J1 + I].conjugate();
          A[DPOS - NB - I][J1 + I] = Complex.zero;
        } // 30
        CTMP.value = A[DPOS - NB][J1];
        zlarfg(LM, CTMP, V(VPOS + 1), 1, TAU(TAUPOS));
        A[DPOS - NB][J1] = CTMP.value;

        zlarfx('Right', LN - 1, LM, V(VPOS), TAU[TAUPOS], A(DPOS - NB + 1, J1),
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

      V[VPOS] = Complex.one;
      for (I = 1; I <= LM - 1; I++) {
        // 20
        V[VPOS + I] = A[OFDPOS + I][ST - 1];
        A[OFDPOS + I][ST - 1] = Complex.zero;
      } // 20
      zlarfg(LM, A(OFDPOS, ST - 1), V(VPOS + 1), 1, TAU(TAUPOS));

      LM = ED - ST + 1;

      zlarfy(UPLO, LM, V(VPOS), 1, TAU[TAUPOS].conjugate(), A(DPOS, ST),
          LDA - 1, WORK);
    }

    if (TTYPE == 3) {
      LM = ED - ST + 1;

      zlarfy(UPLO, LM, V(VPOS), 1, TAU[TAUPOS].conjugate(), A(DPOS, ST),
          LDA - 1, WORK);
    }

    if (TTYPE == 2) {
      J1 = ED + 1;
      J2 = min(ED + NB, N);
      LN = ED - ST + 1;
      LM = J2 - J1 + 1;

      if (LM > 0) {
        zlarfx('Right', LM, LN, V(VPOS), TAU[TAUPOS], A(DPOS + NB, ST), LDA - 1,
            WORK);

        if (WANTZ) {
          VPOS = ((SWEEP - 1) % 2) * N + J1;
          TAUPOS = ((SWEEP - 1) % 2) * N + J1;
        } else {
          VPOS = ((SWEEP - 1) % 2) * N + J1;
          TAUPOS = ((SWEEP - 1) % 2) * N + J1;
        }

        V[VPOS] = Complex.one;
        for (I = 1; I <= LM - 1; I++) {
          // 40
          V[VPOS + I] = A[DPOS + NB + I][ST];
          A[DPOS + NB + I][ST] = Complex.zero;
        } // 40
        zlarfg(LM, A(DPOS + NB, ST), V(VPOS + 1), 1, TAU(TAUPOS));

        zlarfx('Left', LM, LN - 1, V(VPOS), TAU[TAUPOS].conjugate(),
            A(DPOS + NB - 1, ST + 1), LDA - 1, WORK);
      }
    }
  }
}
