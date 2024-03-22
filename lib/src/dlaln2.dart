import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dladiv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dlaln2(
  final bool LTRANS,
  final int NA,
  final int NW,
  final double SMIN,
  final double CA,
  final Matrix<double> A_,
  final int LDA,
  final double D1,
  final double D2,
  final Matrix<double> B_,
  final int LDB,
  final double WR,
  final double WI,
  final Matrix<double> X_,
  final int LDX,
  final Box<double> SCALE,
  final Box<double> XNORM,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  const ZERO = 0.0, ONE = 1.0;
  const TWO = 2.0;
  int ICMAX, J;
  double BBND,
      BI1,
      BI2,
      BIGNUM,
      BNORM,
      BR1,
      BR2,
      CI21,
      CI22,
      CMAX,
      CNORM,
      CR21,
      CR22,
      CSI,
      CSR,
      LI21,
      LR21,
      SMINI,
      SMLNUM,
      TEMP,
      U22ABS,
      UI11,
      UI11R,
      UI12,
      UI12S,
      UI22,
      UR11,
      UR11R,
      UR12,
      UR12S,
      UR22,
      XI1 = 0,
      XR1 = 0;
  final CI = Matrix<double>(2, 2), CR = Matrix<double>(2, 2);
  final CIV = CI(1, 1).asArray(), CRV = CR(1, 1).asArray();
  final ZSWAP = Array.fromList([false, false, true, true]);
  final RSWAP = Array.fromList([false, true, false, true]);
  final IPIVOT = Matrix<int>.fromList([
    [1, 2, 3, 4],
    [2, 1, 4, 3],
    [3, 4, 1, 2],
    [4, 3, 2, 1],
  ]);
  final XI2 = Box(0.0), XR2 = Box(0.0);

  // Compute BIGNUM

  SMLNUM = TWO * dlamch('Safe minimum');
  BIGNUM = ONE / SMLNUM;
  SMINI = max(SMIN, SMLNUM);

  // Don't check for input errors

  INFO.value = 0;

  // Standard Initializations

  SCALE.value = ONE;

  if (NA == 1) {
    // 1 x 1  (i.e., scalar) system   C X = B

    if (NW == 1) {
      // Real 1x1 system.

      // C = ca A - w D

      CSR = CA * A[1][1] - WR * D1;
      CNORM = CSR.abs();

      // If | C | < SMINI, use C = SMINI

      if (CNORM < SMINI) {
        CSR = SMINI;
        CNORM = SMINI;
        INFO.value = 1;
      }

      // Check scaling for  X = B / C

      BNORM = B[1][1].abs();
      if (CNORM < ONE && BNORM > ONE) {
        if (BNORM > BIGNUM * CNORM) SCALE.value = ONE / BNORM;
      }

      // Compute X

      X[1][1] = (B[1][1] * SCALE.value) / CSR;
      XNORM.value = X[1][1].abs();
    } else {
      // Complex 1x1 system (w is complex)

      // C = ca A - w D

      CSR = CA * A[1][1] - WR * D1;
      CSI = -WI * D1;
      CNORM = CSR.abs() + CSI.abs();

      // If | C | < SMINI, use C = SMINI

      if (CNORM < SMINI) {
        CSR = SMINI;
        CSI = ZERO;
        CNORM = SMINI;
        INFO.value = 1;
      }

      // Check scaling for  X = B / C

      BNORM = B[1][1].abs() + B[1][2].abs();
      if (CNORM < ONE && BNORM > ONE) {
        if (BNORM > BIGNUM * CNORM) SCALE.value = ONE / BNORM;
      }

      // Compute X

      dladiv(SCALE.value * B[1][1], SCALE.value * B[1][2], CSR, CSI,
          X.box(1, 1), X.box(1, 2));
      XNORM.value = X[1][1].abs() + X[1][2].abs();
    }
  } else {
    // 2x2 System

    // Compute the real part of  C = ca A - w D  (or  ca A**T - w D )

    CR[1][1] = CA * A[1][1] - WR * D1;
    CR[2][2] = CA * A[2][2] - WR * D2;
    if (LTRANS) {
      CR[1][2] = CA * A[2][1];
      CR[2][1] = CA * A[1][2];
    } else {
      CR[2][1] = CA * A[2][1];
      CR[1][2] = CA * A[1][2];
    }

    if (NW == 1) {
      // Real 2x2 system  (w is real)

      // Find the largest element in C

      CMAX = ZERO;
      ICMAX = 0;

      for (J = 1; J <= 4; J++) {
        if (CRV[J].abs() > CMAX) {
          CMAX = CRV[J].abs();
          ICMAX = J;
        }
      }

      // If norm(C) < SMINI, use SMINI*identity.

      if (CMAX < SMINI) {
        BNORM = max(B[1][1].abs(), B[2][1].abs());
        if (SMINI < ONE && BNORM > ONE) {
          if (BNORM > BIGNUM * SMINI) SCALE.value = ONE / BNORM;
        }
        TEMP = SCALE.value / SMINI;
        X[1][1] = TEMP * B[1][1];
        X[2][1] = TEMP * B[2][1];
        XNORM.value = TEMP * BNORM;
        INFO.value = 1;
        return;
      }

      // Gaussian elimination with complete pivoting.

      UR11 = CRV[ICMAX];
      CR21 = CRV[IPIVOT[2][ICMAX]];
      UR12 = CRV[IPIVOT[3][ICMAX]];
      CR22 = CRV[IPIVOT[4][ICMAX]];
      UR11R = ONE / UR11;
      LR21 = UR11R * CR21;
      UR22 = CR22 - UR12 * LR21;

      // If smaller pivot < SMINI, use SMINI

      if (UR22.abs() < SMINI) {
        UR22 = SMINI;
        INFO.value = 1;
      }
      if (RSWAP[ICMAX]) {
        BR1 = B[2][1];
        BR2 = B[1][1];
      } else {
        BR1 = B[1][1];
        BR2 = B[2][1];
      }
      BR2 -= LR21 * BR1;
      BBND = max((BR1 * (UR22 * UR11R).abs()), BR2.abs());
      if (BBND > ONE && UR22.abs() < ONE) {
        if (BBND >= BIGNUM * UR22.abs()) SCALE.value = ONE / BBND;
      }

      XR2.value = (BR2 * SCALE.value) / UR22;
      XR1 = (SCALE.value * BR1) * UR11R - XR2.value * (UR11R * UR12);
      if (ZSWAP[ICMAX]) {
        X[1][1] = XR2.value;
        X[2][1] = XR1;
      } else {
        X[1][1] = XR1;
        X[2][1] = XR2.value;
      }
      XNORM.value = max(XR1.abs(), XR2.value.abs());

      // Further scaling if  norm(A) norm(X) > overflow

      if (XNORM.value > ONE && CMAX > ONE) {
        if (XNORM.value > BIGNUM / CMAX) {
          TEMP = CMAX / BIGNUM;
          X[1][1] = TEMP * X[1][1];
          X[2][1] = TEMP * X[2][1];
          XNORM.value = TEMP * XNORM.value;
          SCALE.value = TEMP * SCALE.value;
        }
      }
    } else {
      // Complex 2x2 system  (w is complex)

      // Find the largest element in C

      CI[1][1] = -WI * D1;
      CI[2][1] = ZERO;
      CI[1][2] = ZERO;
      CI[2][2] = -WI * D2;
      CMAX = ZERO;
      ICMAX = 0;

      for (J = 1; J <= 4; J++) {
        if (CRV[J].abs() + CIV[J].abs() > CMAX) {
          CMAX = CRV[J].abs() + CIV[J].abs();
          ICMAX = J;
        }
      }

      // If norm(C) < SMINI, use SMINI*identity.

      if (CMAX < SMINI) {
        BNORM =
            max(B[1][1].abs() + B[1][2].abs(), B[2][1].abs() + B[2][2]).abs();
        if (SMINI < ONE && BNORM > ONE) {
          if (BNORM > BIGNUM * SMINI) SCALE.value = ONE / BNORM;
        }
        TEMP = SCALE.value / SMINI;
        X[1][1] = TEMP * B[1][1];
        X[2][1] = TEMP * B[2][1];
        X[1][2] = TEMP * B[1][2];
        X[2][2] = TEMP * B[2][2];
        XNORM.value = TEMP * BNORM;
        INFO.value = 1;
        return;
      }

      // Gaussian elimination with complete pivoting.

      UR11 = CRV[ICMAX];
      UI11 = CIV[ICMAX];
      CR21 = CRV[IPIVOT[2][ICMAX]];
      CI21 = CIV[IPIVOT[2][ICMAX]];
      UR12 = CRV[IPIVOT[3][ICMAX]];
      UI12 = CIV[IPIVOT[3][ICMAX]];
      CR22 = CRV[IPIVOT[4][ICMAX]];
      CI22 = CIV[IPIVOT[4][ICMAX]];
      if (ICMAX == 1 || ICMAX == 4) {
        // Code when off-diagonals of pivoted C are real

        if (UR11.abs() > UI11.abs()) {
          TEMP = UI11 / UR11;
          UR11R = ONE / (UR11 * (ONE + pow(TEMP, 2)));
          UI11R = -TEMP * UR11R;
        } else {
          TEMP = UR11 / UI11;
          UI11R = -ONE / (UI11 * (ONE + pow(TEMP, 2)));
          UR11R = -TEMP * UI11R;
        }
        LR21 = CR21 * UR11R;
        LI21 = CR21 * UI11R;
        UR12S = UR12 * UR11R;
        UI12S = UR12 * UI11R;
        UR22 = CR22 - UR12 * LR21;
        UI22 = CI22 - UR12 * LI21;
      } else {
        // Code when diagonals of pivoted C are real

        UR11R = ONE / UR11;
        UI11R = ZERO;
        LR21 = CR21 * UR11R;
        LI21 = CI21 * UR11R;
        UR12S = UR12 * UR11R;
        UI12S = UI12 * UR11R;
        UR22 = CR22 - UR12 * LR21 + UI12 * LI21;
        UI22 = -UR12 * LI21 - UI12 * LR21;
      }
      U22ABS = UR22.abs() + UI22.abs();

      // If smaller pivot < SMINI, use SMINI

      if (U22ABS < SMINI) {
        UR22 = SMINI;
        UI22 = ZERO;
        INFO.value = 1;
      }
      if (RSWAP[ICMAX]) {
        BR2 = B[1][1];
        BR1 = B[2][1];
        BI2 = B[1][2];
        BI1 = B[2][2];
      } else {
        BR1 = B[1][1];
        BR2 = B[2][1];
        BI1 = B[1][2];
        BI2 = B[2][2];
      }
      BR2 -= LR21 * BR1 - LI21 * BI1;
      BI2 -= LI21 * BR1 + LR21 * BI1;
      BBND = max(
          (BR1.abs() + BI1.abs()) * (U22ABS * (UR11R.abs() + UI11R.abs())),
          BR2.abs() + BI2.abs());
      if (BBND > ONE && U22ABS < ONE) {
        if (BBND >= BIGNUM * U22ABS) {
          SCALE.value = ONE / BBND;
          BR1 = SCALE.value * BR1;
          BI1 = SCALE.value * BI1;
          BR2 = SCALE.value * BR2;
          BI2 = SCALE.value * BI2;
        }
      }

      dladiv(BR2, BI2, UR22, UI22, XR2, XI2);
      XR1 = UR11R * BR1 - UI11R * BI1 - UR12S * XR2.value + UI12S * XI2.value;
      XI1 = UI11R * BR1 + UR11R * BI1 - UI12S * XR2.value - UR12S * XI2.value;
      if (ZSWAP[ICMAX]) {
        X[1][1] = XR2.value;
        X[2][1] = XR1;
        X[1][2] = XI2.value;
        X[2][2] = XI1;
      } else {
        X[1][1] = XR1;
        X[2][1] = XR2.value;
        X[1][2] = XI1;
        X[2][2] = XI2.value;
      }
      XNORM.value =
          max(XR1.abs() + XI1.abs(), XR2.value.abs() + XI2.value.abs());

      // Further scaling if  norm(A) norm(X) > overflow

      if (XNORM.value > ONE && CMAX > ONE) {
        if (XNORM.value > BIGNUM / CMAX) {
          TEMP = CMAX / BIGNUM;
          X[1][1] = TEMP * X[1][1];
          X[2][1] = TEMP * X[2][1];
          X[1][2] = TEMP * X[1][2];
          X[2][2] = TEMP * X[2][2];
          XNORM.value = TEMP * XNORM.value;
          SCALE.value = TEMP * SCALE.value;
        }
      }
    }
  }
}
