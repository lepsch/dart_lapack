import 'dart:math';

import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zgehrd.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlahqr.dart';
import 'package:lapack/src/zlarf.dart';
import 'package:lapack/src/zlarfg.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/ztrexc.dart';
import 'package:lapack/src/zunmhr.dart';

void zlaqr2(
  final bool WANTT,
  final bool WANTZ,
  final int N,
  final int KTOP,
  final int KBOT,
  final int NW,
  final Matrix<Complex> H_,
  final int LDH,
  final int ILOZ,
  final int IHIZ,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Box<int> NS,
  final Box<int> ND,
  final Array<Complex> SH_,
  final Matrix<Complex> V_,
  final int LDV,
  final int NH,
  final Matrix<Complex> T_,
  final int LDT,
  final int NV,
  final Matrix<Complex> WV_,
  final int LDWV,
  final Array<Complex> WORK_,
  final int LWORK,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final H = H_.having(ld: LDH);
  final Z = Z_.having(ld: LDZ);
  final V = V_.having(ld: LDV);
  final T = T_.having(ld: LDT);
  final WV = WV_.having(ld: LDWV);
  final SH = SH_.having();
  final WORK = WORK_.having();
  const RZERO = 0.0;
  Complex S;
  double FOO,
      // SAFMAX,
      SAFMIN,
      SMLNUM,
      ULP;
  int I,
      IFST,
      ILST,
      J,
      JW,
      KCOL,
      KLN,
      KNT,
      KROW,
      KWTOP,
      LTOP,
      LWK1,
      LWK2,
      LWKOPT;
  final INFO = Box(0), INFQR = Box(0);
  final BETA = Box(Complex.zero), TAU = Box(Complex.zero);

  double CABS1(Complex CDUM) => CDUM.toDouble().abs() + CDUM.imaginary.abs();

  // ==== Estimate optimal workspace. ====

  JW = min(NW, KBOT - KTOP + 1);
  if (JW <= 2) {
    LWKOPT = 1;
  } else {
    // ==== Workspace query call to ZGEHRD ====

    zgehrd(JW, 1, JW - 1, T, LDT, WORK, WORK, -1, INFO);
    LWK1 = WORK[1].toInt();

    // ==== Workspace query call to ZUNMHR ====

    zunmhr('R', 'N', JW, JW, 1, JW - 1, T, LDT, WORK, V, LDV, WORK, -1, INFO);
    LWK2 = WORK[1].toInt();

    // ==== Optimal workspace ====

    LWKOPT = JW + max(LWK1, LWK2);
  }

  // ==== Quick return in case of workspace query. ====

  if (LWORK == -1) {
    WORK[1] = LWKOPT.toComplex();
    return;
  }

  // ==== Nothing to do ...
  // ... for an empty active block ... ====
  NS.value = 0;
  ND.value = 0;
  WORK[1] = Complex.one;
  if (KTOP > KBOT) return;
  // ... nor for an empty deflation window. ====
  if (NW < 1) return;

  // ==== Machine constants ====

  SAFMIN = dlamch('SAFE MINIMUM');
  // SAFMAX = RONE / SAFMIN;
  ULP = dlamch('PRECISION');
  SMLNUM = SAFMIN * (N.toDouble() / ULP);

  // ==== Setup deflation window ====

  JW = min(NW, KBOT - KTOP + 1);
  KWTOP = KBOT - JW + 1;
  if (KWTOP == KTOP) {
    S = Complex.zero;
  } else {
    S = H[KWTOP][KWTOP - 1];
  }

  if (KBOT == KWTOP) {
    // ==== 1-by-1 deflation window: not much to do ====

    SH[KWTOP] = H[KWTOP][KWTOP];
    NS.value = 1;
    ND.value = 0;
    if (CABS1(S) <= max(SMLNUM, ULP * CABS1(H[KWTOP][KWTOP]))) {
      NS.value = 0;
      ND.value = 1;
      if (KWTOP > KTOP) H[KWTOP][KWTOP - 1] = Complex.zero;
    }
    WORK[1] = Complex.one;
    return;
  }

  // ==== Convert to spike-triangular form.  (In case of a
  // .    rare QR failure, this routine continues to do
  // .    aggressive early deflation using that part of
  // .    the deflation window that converged using INFQR
  // .    here and there to keep track.) ====

  zlacpy('U', JW, JW, H(KWTOP, KWTOP), LDH, T, LDT);
  zcopy(JW - 1, H(KWTOP + 1, KWTOP).asArray(), LDH + 1, T(2, 1).asArray(),
      LDT + 1);

  zlaset('A', JW, JW, Complex.zero, Complex.one, V, LDV);
  zlahqr(true, true, JW, 1, JW, T, LDT, SH(KWTOP), 1, JW, V, LDV, INFQR);

  // ==== Deflation detection loop ====

  NS.value = JW;
  ILST = INFQR.value + 1;
  for (KNT = INFQR.value + 1; KNT <= JW; KNT++) {
    // 10

    // ==== Small spike tip deflation test ====

    FOO = CABS1(T[NS.value][NS.value]);
    if (FOO == RZERO) FOO = CABS1(S);
    if (CABS1(S) * CABS1(V[1][NS.value]) <= max(SMLNUM, ULP * FOO)) {
      // ==== One more converged eigenvalue ====

      NS.value = NS.value - 1;
    } else {
      // ==== One undeflatable eigenvalue.  Move it up out of the
      // .    way.   (ZTREXC can not fail in this case.) ====

      IFST = NS.value;
      ztrexc('V', JW, T, LDT, V, LDV, IFST, ILST, INFO);
      ILST++;
    }
  } // 10

  // ==== Return to Hessenberg form ====

  if (NS.value == 0) S = Complex.zero;

  if (NS.value < JW) {
    // ==== sorting the diagonal of T improves accuracy for
    // .    graded matrices.  ====

    for (I = INFQR.value + 1; I <= NS.value; I++) {
      // 30
      IFST = I;
      for (J = I + 1; J <= NS.value; J++) {
        // 20
        if (CABS1(T[J][J]) > CABS1(T[IFST][IFST])) IFST = J;
      } // 20
      ILST = I;
      if (IFST != ILST) ztrexc('V', JW, T, LDT, V, LDV, IFST, ILST, INFO);
    } // 30
  }

  // ==== Restore shift/eigenvalue array from T ====

  for (I = INFQR.value + 1; I <= JW; I++) {
    // 40
    SH[KWTOP + I - 1] = T[I][I];
  } // 40

  if (NS.value < JW || S == Complex.zero) {
    if (NS.value > 1 && S != Complex.zero) {
      // ==== Reflect spike back into lower triangle ====

      zcopy(NS.value, V.asArray(), LDV, WORK, 1);
      for (I = 1; I <= NS.value; I++) {
        // 50
        WORK[I] = WORK[I].conjugate();
      } // 50
      BETA.value = WORK[1];
      zlarfg(NS.value, BETA, WORK(2), 1, TAU);
      WORK[1] = Complex.one;

      zlaset('L', JW - 2, JW - 2, Complex.zero, Complex.zero, T(3, 1), LDT);

      zlarf('L', NS.value, JW, WORK, 1, TAU.value.conjugate(), T, LDT,
          WORK(JW + 1));
      zlarf('R', NS.value, NS.value, WORK, 1, TAU.value, T, LDT, WORK(JW + 1));
      zlarf('R', JW, NS.value, WORK, 1, TAU.value, V, LDV, WORK(JW + 1));

      zgehrd(JW, 1, NS.value, T, LDT, WORK, WORK(JW + 1), LWORK - JW, INFO);
    }

    // ==== Copy updated reduced window into place ====

    if (KWTOP > 1) H[KWTOP][KWTOP - 1] = S * V[1][1].conjugate();
    zlacpy('U', JW, JW, T, LDT, H(KWTOP, KWTOP), LDH);
    zcopy(JW - 1, T(2, 1).asArray(), LDT + 1, H(KWTOP + 1, KWTOP).asArray(),
        LDH + 1);

    // ==== Accumulate orthogonal matrix in order update
    // .    H and Z, if requested.  ====

    if (NS.value > 1 && S != Complex.zero) {
      zunmhr('R', 'N', JW, NS.value, 1, NS.value, T, LDT, WORK, V, LDV,
          WORK(JW + 1), LWORK - JW, INFO);
    }

    // ==== Update vertical slab in H ====

    if (WANTT) {
      LTOP = 1;
    } else {
      LTOP = KTOP;
    }
    for (KROW = LTOP;
        NV < 0 ? KROW >= KWTOP - 1 : KROW <= KWTOP - 1;
        KROW += NV) {
      // 60
      KLN = min(NV, KWTOP - KROW);
      zgemm('N', 'N', KLN, JW, JW, Complex.one, H(KROW, KWTOP), LDH, V, LDV,
          Complex.zero, WV, LDWV);
      zlacpy('A', KLN, JW, WV, LDWV, H(KROW, KWTOP), LDH);
    } // 60

    // ==== Update horizontal slab in H ====

    if (WANTT) {
      for (KCOL = KBOT + 1; NH < 0 ? KCOL >= N : KCOL <= N; KCOL += NH) {
        // 70
        KLN = min(NH, N - KCOL + 1);
        zgemm('C', 'N', JW, KLN, JW, Complex.one, V, LDV, H(KWTOP, KCOL), LDH,
            Complex.zero, T, LDT);
        zlacpy('A', JW, KLN, T, LDT, H(KWTOP, KCOL), LDH);
      } // 70
    }

    // ==== Update vertical slab in Z ====

    if (WANTZ) {
      for (KROW = ILOZ; NV < 0 ? KROW >= IHIZ : KROW <= IHIZ; KROW += NV) {
        // 80
        KLN = min(NV, IHIZ - KROW + 1);
        zgemm('N', 'N', KLN, JW, JW, Complex.one, Z(KROW, KWTOP), LDZ, V, LDV,
            Complex.zero, WV, LDWV);
        zlacpy('A', KLN, JW, WV, LDWV, Z(KROW, KWTOP), LDZ);
      } // 80
    }
  }

  // ==== Return the number of deflations ... ====

  ND.value = JW - NS.value;

  // ==== ... and the number of shifts. (Subtracting
  // .    INFQR.value from the spike length takes care
  // .    of the case of a rare QR failure while
  // .    calculating eigenvalues of the deflation
  // .    window.)  ====

  NS.value = NS.value - INFQR.value;

  // ==== Return optimal workspace. ====

  WORK[1] = LWKOPT.toComplex();
}
