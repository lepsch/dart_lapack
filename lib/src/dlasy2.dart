import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dlasy2(
  final bool LTRANL,
  final bool LTRANR,
  final int ISGN,
  final int N1,
  final int N2,
  final Matrix<double> TL_,
  final int LDTL,
  final Matrix<double> TR_,
  final int LDTR,
  final Matrix<double> B_,
  final int LDB,
  final Box<double> SCALE,
  final Matrix<double> X_,
  final int LDX,
  final Box<double> XNORM,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final TL = TL_.having(ld: LDTL);
  final TR = TR_.having(ld: LDTR);
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  const ZERO = 0.0, ONE = 1.0;
  const TWO = 2.0, HALF = 0.5, EIGHT = 8.0;
  bool BSWAP, XSWAP;
  int I, IP, IPIV, IPSV = 0, J, JP, JPSV = 0, K;
  double BET,
      EPS,
      GAM,
      L21,
      SGN,
      SMIN = 0,
      SMLNUM,
      TAU1,
      TEMP,
      U11,
      U12,
      U22,
      XMAX;
  final JPIV = Array<int>(4);
  final BTMP = Array<double>(4), TMP = Array<double>(4), X2 = Array<double>(2);
  final T16 = Matrix<double>(4, 4);
  final LOCU12 = Array.fromList([3, 4, 1, 2]);
  final LOCL21 = Array.fromList([2, 1, 4, 3]);
  final LOCU22 = Array.fromList([4, 3, 2, 1]);
  final XSWPIV = Array.fromList([false, false, true, true]);
  final BSWPIV = Array.fromList([false, true, false, true]);

  // Do not check the input parameters for errors

  INFO.value = 0;

  // Quick return if possible

  if (N1 == 0 || N2 == 0) return;

  // Set constants to control overflow

  EPS = dlamch('P');
  SMLNUM = dlamch('S') / EPS;
  SGN = ISGN.toDouble();

  K = N1 + N1 + N2 - 2;
  // GO TO ( 10, 20, 30, 50 )K;
  if (K == 1) {
    // 1 by 1: TL11*X + SGN*X*TR11 = B11

    TAU1 = TL[1][1] + SGN * TR[1][1];
    BET = TAU1.abs();
    if (BET <= SMLNUM) {
      TAU1 = SMLNUM;
      BET = SMLNUM;
      INFO.value = 1;
    }

    SCALE.value = ONE;
    GAM = B[1][1].abs();
    if (SMLNUM * GAM > BET) SCALE.value = ONE / GAM;

    X[1][1] = (B[1][1] * SCALE.value) / TAU1;
    XNORM.value = X[1][1].abs();
    return;
  }

  if (K == 2) {
    // 1 by 2:
    // TL11*[X11 X12] + ISGN*[X11 X12]*op[TR11 TR12]  = [B11 B12]
    //                                   [TR21 TR22]

    SMIN = max(
        EPS *
            max(
              TL[1][1].abs(),
              max(
                TR[1][1].abs(),
                max(
                  TR[1][2].abs(),
                  max(TR[2][1].abs(), TR[2][2]),
                ),
              ),
            ).abs(),
        SMLNUM);
    TMP[1] = TL[1][1] + SGN * TR[1][1];
    TMP[4] = TL[1][1] + SGN * TR[2][2];
    if (LTRANR) {
      TMP[2] = SGN * TR[2][1];
      TMP[3] = SGN * TR[1][2];
    } else {
      TMP[2] = SGN * TR[1][2];
      TMP[3] = SGN * TR[2][1];
    }
    BTMP[1] = B[1][1];
    BTMP[2] = B[1][2];
  } else if (K == 3) {
    // 2 by 1:
    //      op[TL11 TL12]*[X11] + ISGN* [X11]*TR11  = [B11]
    //        [TL21 TL22] [X21]         [X21]         [B21]

    SMIN = max(
        EPS *
            max(
              TR[1][1].abs(),
              max(
                TL[1][1].abs(),
                max(
                  TL[1][2].abs(),
                  max(TL[2][1].abs(), TL[2][2]),
                ),
              ),
            ).abs(),
        SMLNUM);
    TMP[1] = TL[1][1] + SGN * TR[1][1];
    TMP[4] = TL[2][2] + SGN * TR[1][1];
    if (LTRANL) {
      TMP[2] = TL[1][2];
      TMP[3] = TL[2][1];
    } else {
      TMP[2] = TL[2][1];
      TMP[3] = TL[1][2];
    }
    BTMP[1] = B[1][1];
    BTMP[2] = B[2][1];
  }

  if (K != 4) {
    // Solve 2 by 2 system using complete pivoting.
    // Set pivots less than SMIN to SMIN.

    IPIV = idamax(4, TMP, 1);
    U11 = TMP[IPIV];
    if (U11.abs() <= SMIN) {
      INFO.value = 1;
      U11 = SMIN;
    }
    U12 = TMP[LOCU12[IPIV]];
    L21 = TMP[LOCL21[IPIV]] / U11;
    U22 = TMP[LOCU22[IPIV]] - U12 * L21;
    XSWAP = XSWPIV[IPIV];
    BSWAP = BSWPIV[IPIV];
    if (U22.abs() <= SMIN) {
      INFO.value = 1;
      U22 = SMIN;
    }
    if (BSWAP) {
      TEMP = BTMP[2];
      BTMP[2] = BTMP[1] - L21 * TEMP;
      BTMP[1] = TEMP;
    } else {
      BTMP[2] -= L21 * BTMP[1];
    }
    SCALE.value = ONE;
    if ((TWO * SMLNUM) * BTMP[2].abs() > U22.abs() ||
        (TWO * SMLNUM) * BTMP[1].abs() > U11.abs()) {
      SCALE.value = HALF / max(BTMP[1].abs(), BTMP[2]).abs();
      BTMP[1] *= SCALE.value;
      BTMP[2] *= SCALE.value;
    }
    X2[2] = BTMP[2] / U22;
    X2[1] = BTMP[1] / U11 - (U12 / U11) * X2[2];
    if (XSWAP) {
      TEMP = X2[2];
      X2[2] = X2[1];
      X2[1] = TEMP;
    }
    X[1][1] = X2[1];
    if (N1 == 1) {
      X[1][2] = X2[2];
      XNORM.value = X[1][1].abs() + X[1][2].abs();
    } else {
      X[2][1] = X2[2];
      XNORM.value = max(X[1][1].abs(), X[2][1]).abs();
    }
    return;
  }

  // 2 by 2:
  // op[TL11 TL12]*[X11 X12] +ISGN* [X11 X12]*op[TR11 TR12] = [B11 B12]
  //   [TL21 TL22] [X21 X22]        [X21 X22]   [TR21 TR22]   [B21 B22]

  // Solve equivalent 4 by 4 system using complete pivoting.
  // Set pivots less than SMIN to SMIN.

  SMIN = max(
      TR[1][1].abs(), max(TR[1][2].abs(), max(TR[2][1].abs(), TR[2][2].abs())));
  SMIN = max(
      SMIN,
      max(
        TL[1][1].abs(),
        max(TL[1][2].abs(), max(TL[2][1].abs(), TL[2][2].abs())),
      ));
  SMIN = max(EPS * SMIN, SMLNUM);
  BTMP[1] = ZERO;
  dcopy(16, BTMP, 0, T16.asArray(), 1);
  T16[1][1] = TL[1][1] + SGN * TR[1][1];
  T16[2][2] = TL[2][2] + SGN * TR[1][1];
  T16[3][3] = TL[1][1] + SGN * TR[2][2];
  T16[4][4] = TL[2][2] + SGN * TR[2][2];
  if (LTRANL) {
    T16[1][2] = TL[2][1];
    T16[2][1] = TL[1][2];
    T16[3][4] = TL[2][1];
    T16[4][3] = TL[1][2];
  } else {
    T16[1][2] = TL[1][2];
    T16[2][1] = TL[2][1];
    T16[3][4] = TL[1][2];
    T16[4][3] = TL[2][1];
  }
  if (LTRANR) {
    T16[1][3] = SGN * TR[1][2];
    T16[2][4] = SGN * TR[1][2];
    T16[3][1] = SGN * TR[2][1];
    T16[4][2] = SGN * TR[2][1];
  } else {
    T16[1][3] = SGN * TR[2][1];
    T16[2][4] = SGN * TR[2][1];
    T16[3][1] = SGN * TR[1][2];
    T16[4][2] = SGN * TR[1][2];
  }
  BTMP[1] = B[1][1];
  BTMP[2] = B[2][1];
  BTMP[3] = B[1][2];
  BTMP[4] = B[2][2];

  // Perform elimination

  for (I = 1; I <= 3; I++) {
    XMAX = ZERO;
    for (IP = I; IP <= 4; IP++) {
      for (JP = I; JP <= 4; JP++) {
        if ((T16[IP][JP].abs()) >= XMAX) {
          XMAX = (T16[IP][JP].abs());
          IPSV = IP;
          JPSV = JP;
        }
      }
    }
    if (IPSV != I) {
      dswap(4, T16(IPSV, 1).asArray(), 4, T16(I, 1).asArray(), 4);
      TEMP = BTMP[I];
      BTMP[I] = BTMP[IPSV];
      BTMP[IPSV] = TEMP;
    }
    if (JPSV != I) dswap(4, T16(1, JPSV).asArray(), 1, T16(1, I).asArray(), 1);
    JPIV[I] = JPSV;
    if (T16[I][I].abs() < SMIN) {
      INFO.value = 1;
      T16[I][I] = SMIN;
    }
    for (J = I + 1; J <= 4; J++) {
      T16[J][I] /= T16[I][I];
      BTMP[J] -= T16[J][I] * BTMP[I];
      for (K = I + 1; K <= 4; K++) {
        T16[J][K] -= T16[J][I] * T16[I][K];
      }
    }
  }
  if (T16[4][4].abs() < SMIN) {
    INFO.value = 1;
    T16[4][4] = SMIN;
  }
  SCALE.value = ONE;
  if ((EIGHT * SMLNUM) * BTMP[1].abs() > T16[1][1].abs() ||
      (EIGHT * SMLNUM) * BTMP[2].abs() > T16[2][2].abs() ||
      (EIGHT * SMLNUM) * BTMP[3].abs() > T16[3][3].abs() ||
      (EIGHT * SMLNUM) * BTMP[4].abs() > T16[4][4].abs()) {
    SCALE.value = (ONE / EIGHT) /
        max(BTMP[1].abs(),
            max(BTMP[2].abs(), max(BTMP[3].abs(), BTMP[4].abs())));
    BTMP[1] *= SCALE.value;
    BTMP[2] *= SCALE.value;
    BTMP[3] *= SCALE.value;
    BTMP[4] *= SCALE.value;
  }
  for (I = 1; I <= 4; I++) {
    K = 5 - I;
    TEMP = ONE / T16[K][K];
    TMP[K] = BTMP[K] * TEMP;
    for (J = K + 1; J <= 4; J++) {
      TMP[K] -= (TEMP * T16[K][J]) * TMP[J];
    }
  }
  for (I = 1; I <= 3; I++) {
    if (JPIV[4 - I] != 4 - I) {
      TEMP = TMP[4 - I];
      TMP[4 - I] = TMP[JPIV[4 - I]];
      TMP[JPIV[4 - I]] = TEMP;
    }
  }
  X[1][1] = TMP[1];
  X[2][1] = TMP[2];
  X[1][2] = TMP[3];
  X[2][2] = TMP[4];
  XNORM.value = max(TMP[1].abs() + TMP[3].abs(), TMP[2].abs() + TMP[4]).abs();
}
