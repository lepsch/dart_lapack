import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlaln2.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dget31(
  final Box<double> RMAX,
  final Box<int> LMAX,
  final Array<int> NINFO_,
  final Box<int> KNT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final NINFO = NINFO_.having();
  const ZERO = 0.0, HALF = 0.5, ONE = 1.0;
  const TWO = 2.0, THREE = 3.0, FOUR = 4.0;
  const SEVEN = 7.0, TEN = 10.0;
  const TWNONE = 21.0;
  int IA, IB, ICA, ID1, ID2, ISMIN, ITRANS, IWI, IWR, NA, NW;
  double BIGNUM, CA, D1, D2, DEN, EPS, RES, SMIN, SMLNUM, TMP, UNFL, WI, WR;
  final A = Matrix<double>(2, 2),
      B = Matrix<double>(2, 2),
      X = Matrix<double>(2, 2);
  final VAB = Array<double>(3),
      VCA = Array<double>(5),
      VDD = Array<double>(4),
      VSMIN = Array<double>(4),
      VWI = Array<double>(4),
      VWR = Array<double>(4);
  final LTRANS = Array.fromList([false, true], offset: zeroIndexedArrayOffset);
  final SCALE = Box(0.0), XNORM = Box(0.0);
  final INFO = Box(0);

  // Get machine parameters

  EPS = dlamch('P');
  UNFL = dlamch('U');
  SMLNUM = dlamch('S') / EPS;
  BIGNUM = ONE / SMLNUM;

  // Set up test case parameters

  VSMIN[1] = SMLNUM;
  VSMIN[2] = EPS;
  VSMIN[3] = ONE / (TEN * TEN);
  VSMIN[4] = ONE / EPS;
  VAB[1] = sqrt(SMLNUM);
  VAB[2] = ONE;
  VAB[3] = sqrt(BIGNUM);
  VWR[1] = ZERO;
  VWR[2] = HALF;
  VWR[3] = TWO;
  VWR[4] = ONE;
  VWI[1] = SMLNUM;
  VWI[2] = EPS;
  VWI[3] = ONE;
  VWI[4] = TWO;
  VDD[1] = sqrt(SMLNUM);
  VDD[2] = ONE;
  VDD[3] = TWO;
  VDD[4] = sqrt(BIGNUM);
  VCA[1] = ZERO;
  VCA[2] = sqrt(SMLNUM);
  VCA[3] = EPS;
  VCA[4] = HALF;
  VCA[5] = ONE;

  KNT.value = 0;
  NINFO[1] = 0;
  NINFO[2] = 0;
  LMAX.value = 0;
  RMAX.value = ZERO;

  // Begin test loop

  for (ID1 = 1; ID1 <= 4; ID1++) {
    D1 = VDD[ID1];
    for (ID2 = 1; ID2 <= 4; ID2++) {
      D2 = VDD[ID2];
      for (ICA = 1; ICA <= 5; ICA++) {
        CA = VCA[ICA];
        for (ITRANS = 0; ITRANS <= 1; ITRANS++) {
          for (ISMIN = 1; ISMIN <= 4; ISMIN++) {
            SMIN = VSMIN[ISMIN];

            NA = 1;
            NW = 1;
            for (IA = 1; IA <= 3; IA++) {
              A[1][1] = VAB[IA];
              for (IB = 1; IB <= 3; IB++) {
                B[1][1] = VAB[IB];
                for (IWR = 1; IWR <= 4; IWR++) {
                  if (D1 == ONE && D2 == ONE && CA == ONE) {
                    WR = VWR[IWR] * A[1][1];
                  } else {
                    WR = VWR[IWR];
                  }
                  WI = ZERO;
                  dlaln2(LTRANS[ITRANS], NA, NW, SMIN, CA, A, 2, D1, D2, B, 2,
                      WR, WI, X, 2, SCALE, XNORM, INFO);
                  if (INFO.value < 0) NINFO[1] = NINFO[1] + 1;
                  if (INFO.value > 0) NINFO[2] = NINFO[2] + 1;
                  RES = ((CA * A[1][1] - WR * D1) * X[1][1] -
                          SCALE.value * B[1][1])
                      .abs();
                  if (INFO.value == 0) {
                    DEN = max(EPS * ((CA * A[1][1] - WR * D1) * X[1][1]).abs(),
                        SMLNUM);
                  } else {
                    DEN = max(SMIN * X[1][1].abs(), SMLNUM);
                  }
                  RES = RES / DEN;
                  if (X[1][1].abs() < UNFL &&
                      B[1][1].abs() <=
                          SMLNUM * (CA * A[1][1] - WR * D1).abs()) {
                    RES = ZERO;
                  }
                  if (SCALE.value > ONE) RES = RES + ONE / EPS;
                  RES = RES +
                      (XNORM.value - X[1][1].abs()).abs() /
                          max(SMLNUM, XNORM.value) /
                          EPS;
                  if (INFO.value != 0 && INFO.value != 1) RES = RES + ONE / EPS;
                  KNT.value = KNT.value + 1;
                  if (RES > RMAX.value) {
                    LMAX.value = KNT.value;
                    RMAX.value = RES;
                  }
                }
              }
            }

            NA = 1;
            NW = 2;
            for (IA = 1; IA <= 3; IA++) {
              A[1][1] = VAB[IA];
              for (IB = 1; IB <= 3; IB++) {
                B[1][1] = VAB[IB];
                B[1][2] = -HALF * VAB[IB];
                for (IWR = 1; IWR <= 4; IWR++) {
                  if (D1 == ONE && D2 == ONE && CA == ONE) {
                    WR = VWR[IWR] * A[1][1];
                  } else {
                    WR = VWR[IWR];
                  }
                  for (IWI = 1; IWI <= 4; IWI++) {
                    if (D1 == ONE && D2 == ONE && CA == ONE) {
                      WI = VWI[IWI] * A[1][1];
                    } else {
                      WI = VWI[IWI];
                    }
                    dlaln2(LTRANS[ITRANS], NA, NW, SMIN, CA, A, 2, D1, D2, B, 2,
                        WR, WI, X, 2, SCALE, XNORM, INFO);
                    if (INFO.value < 0) NINFO[1] = NINFO[1] + 1;
                    if (INFO.value > 0) NINFO[2] = NINFO[2] + 1;
                    RES = ((CA * A[1][1] - WR * D1) * X[1][1] +
                            (WI * D1) * X[1][2] -
                            SCALE.value * B[1][1])
                        .abs();
                    RES = RES +
                        ((-WI * D1) * X[1][1] +
                                (CA * A[1][1] - WR * D1) * X[1][2] -
                                SCALE.value * B[1][2])
                            .abs();
                    if (INFO.value == 0) {
                      DEN = max(
                          EPS *
                              (max(
                                    (CA * A[1][1] - WR * D1).abs(),
                                    (D1 * WI).abs(),
                                  ) *
                                  (X[1][1].abs() + X[1][2].abs())),
                          SMLNUM);
                    } else {
                      DEN = max(SMIN * (X[1][1].abs() + X[1][2]).abs(), SMLNUM);
                    }
                    RES = RES / DEN;
                    if (X[1][1].abs() < UNFL &&
                        X[1][2].abs() < UNFL &&
                        B[1][1].abs() <=
                            SMLNUM * (CA * A[1][1] - WR * D1).abs()) RES = ZERO;
                    if (SCALE.value > ONE) RES = RES + ONE / EPS;
                    RES = RES +
                        (XNORM.value - X[1][1].abs() - X[1][2].abs()).abs() /
                            max(SMLNUM, XNORM.value) /
                            EPS;
                    if (INFO.value != 0 && INFO.value != 1) {
                      RES += ONE / EPS;
                    }
                    KNT.value = KNT.value + 1;
                    if (RES > RMAX.value) {
                      LMAX.value = KNT.value;
                      RMAX.value = RES;
                    }
                  }
                }
              }
            }

            NA = 2;
            NW = 1;
            for (IA = 1; IA <= 3; IA++) {
              A[1][1] = VAB[IA];
              A[1][2] = -THREE * VAB[IA];
              A[2][1] = -SEVEN * VAB[IA];
              A[2][2] = TWNONE * VAB[IA];
              for (IB = 1; IB <= 3; IB++) {
                B[1][1] = VAB[IB];
                B[2][1] = -TWO * VAB[IB];
                for (IWR = 1; IWR <= 4; IWR++) {
                  if (D1 == ONE && D2 == ONE && CA == ONE) {
                    WR = VWR[IWR] * A[1][1];
                  } else {
                    WR = VWR[IWR];
                  }
                  WI = ZERO;
                  dlaln2(LTRANS[ITRANS], NA, NW, SMIN, CA, A, 2, D1, D2, B, 2,
                      WR, WI, X, 2, SCALE, XNORM, INFO);
                  if (INFO.value < 0) NINFO[1] = NINFO[1] + 1;
                  if (INFO.value > 0) NINFO[2] = NINFO[2] + 1;
                  if (ITRANS == 1) {
                    TMP = A[1][2];
                    A[1][2] = A[2][1];
                    A[2][1] = TMP;
                  }
                  RES = ((CA * A[1][1] - WR * D1) * X[1][1] +
                          (CA * A[1][2]) * X[2][1] -
                          SCALE.value * B[1][1])
                      .abs();
                  RES = RES +
                      ((CA * A[2][1]) * X[1][1] +
                              (CA * A[2][2] - WR * D2) * X[2][1] -
                              SCALE.value * B[2][1])
                          .abs();
                  if (INFO.value == 0) {
                    DEN = max(
                        EPS *
                            (max(
                                  (CA * A[1][1] - WR * D1).abs() +
                                      (CA * A[1][2]).abs(),
                                  (CA * A[2][1]).abs() +
                                      (CA * A[2][2] - WR * D2).abs(),
                                ) *
                                max(X[1][1].abs(), X[2][1].abs())),
                        SMLNUM);
                  } else {
                    DEN = max(
                        EPS *
                            (max(
                                  SMIN / EPS,
                                  max(
                                    (CA * A[1][1] - WR * D1).abs() +
                                        (CA * A[1][2]).abs(),
                                    (CA * A[2][1]).abs() +
                                        (CA * A[2][2] - WR * D2).abs(),
                                  ),
                                ) *
                                max(X[1][1].abs(), X[2][1].abs())),
                        SMLNUM);
                  }
                  RES = RES / DEN;
                  if (X[1][1].abs() < UNFL &&
                      X[2][1].abs() < UNFL &&
                      B[1][1].abs() + B[2][1].abs() <=
                          SMLNUM *
                              ((CA * A[1][1] - WR * D1).abs() +
                                  (CA * A[1][2]).abs() +
                                  (CA * A[2][1]).abs() +
                                  (CA * A[2][2] - WR * D2).abs())) RES = ZERO;
                  if (SCALE.value > ONE) RES = RES + ONE / EPS;
                  RES = RES +
                      (XNORM.value - max(X[1][1].abs(), X[2][1].abs()).abs()) /
                          max(SMLNUM, XNORM.value) /
                          EPS;
                  if (INFO.value != 0 && INFO.value != 1) RES = RES + ONE / EPS;
                  KNT.value = KNT.value + 1;
                  if (RES > RMAX.value) {
                    LMAX.value = KNT.value;
                    RMAX.value = RES;
                  }
                }
              }
            }

            NA = 2;
            NW = 2;
            for (IA = 1; IA <= 3; IA++) {
              A[1][1] = VAB[IA] * TWO;
              A[1][2] = -THREE * VAB[IA];
              A[2][1] = -SEVEN * VAB[IA];
              A[2][2] = TWNONE * VAB[IA];
              for (IB = 1; IB <= 3; IB++) {
                B[1][1] = VAB[IB];
                B[2][1] = -TWO * VAB[IB];
                B[1][2] = FOUR * VAB[IB];
                B[2][2] = -SEVEN * VAB[IB];
                for (IWR = 1; IWR <= 4; IWR++) {
                  if (D1 == ONE && D2 == ONE && CA == ONE) {
                    WR = VWR[IWR] * A[1][1];
                  } else {
                    WR = VWR[IWR];
                  }
                  for (IWI = 1; IWI <= 4; IWI++) {
                    if (D1 == ONE && D2 == ONE && CA == ONE) {
                      WI = VWI[IWI] * A[1][1];
                    } else {
                      WI = VWI[IWI];
                    }
                    dlaln2(LTRANS[ITRANS], NA, NW, SMIN, CA, A, 2, D1, D2, B, 2,
                        WR, WI, X, 2, SCALE, XNORM, INFO);
                    if (INFO.value < 0) NINFO[1] = NINFO[1] + 1;
                    if (INFO.value > 0) NINFO[2] = NINFO[2] + 1;
                    if (ITRANS == 1) {
                      TMP = A[1][2];
                      A[1][2] = A[2][1];
                      A[2][1] = TMP;
                    }
                    RES = ((CA * A[1][1] - WR * D1) * X[1][1] +
                            (CA * A[1][2]) * X[2][1] +
                            (WI * D1) * X[1][2] -
                            SCALE.value * B[1][1])
                        .abs();
                    RES = RES +
                        ((CA * A[1][1] - WR * D1) * X[1][2] +
                                (CA * A[1][2]) * X[2][2] -
                                (WI * D1) * X[1][1] -
                                SCALE.value * B[1][2])
                            .abs();
                    RES = RES +
                        ((CA * A[2][1]) * X[1][1] +
                                (CA * A[2][2] - WR * D2) * X[2][1] +
                                (WI * D2) * X[2][2] -
                                SCALE.value * B[2][1])
                            .abs();
                    RES = RES +
                        ((CA * A[2][1]) * X[1][2] +
                                (CA * A[2][2] - WR * D2) * X[2][2] -
                                (WI * D2) * X[2][1] -
                                SCALE.value * B[2][2])
                            .abs();
                    if (INFO.value == 0) {
                      DEN = max(
                          EPS *
                              (max(
                                    (CA * A[1][1] - WR * D1).abs() +
                                        (CA * A[1][2]).abs() +
                                        (WI * D1).abs(),
                                    (CA * A[2][1]).abs() +
                                        (CA * A[2][2] - WR * D2).abs() +
                                        (WI * D2).abs(),
                                  ) *
                                  max(
                                    X[1][1].abs() + X[2][1].abs(),
                                    X[1][2].abs() + X[2][2].abs(),
                                  )),
                          SMLNUM);
                    } else {
                      DEN = max(
                          EPS *
                              (max(
                                    SMIN / EPS,
                                    max(
                                      (CA * A[1][1] - WR * D1).abs() +
                                          (CA * A[1][2]).abs() +
                                          (WI * D1).abs(),
                                      (CA * A[2][1]).abs() +
                                          (CA * A[2][2] - WR * D2).abs() +
                                          (WI * D2).abs(),
                                    ),
                                  ) *
                                  max(
                                    X[1][1].abs() + X[2][1].abs(),
                                    X[1][2].abs() + X[2][2].abs(),
                                  )),
                          SMLNUM);
                    }
                    RES = RES / DEN;
                    if (X[1][1].abs() < UNFL &&
                        X[2][1].abs() < UNFL &&
                        X[1][2].abs() < UNFL &&
                        X[2][2].abs() < UNFL &&
                        B[1][1].abs() + B[2][1].abs() <=
                            SMLNUM *
                                ((CA * A[1][1] - WR * D1).abs() +
                                    (CA * A[1][2]).abs() +
                                    (CA * A[2][1]).abs() +
                                    (CA * A[2][2] - WR * D2).abs() +
                                    (WI * D2).abs() +
                                    (WI * D1).abs())) RES = ZERO;
                    if (SCALE.value > ONE) RES = RES + ONE / EPS;
                    RES = RES +
                        (XNORM.value -
                                    max(
                                      X[1][1].abs() + X[1][2].abs(),
                                      X[2][1].abs() + X[2][2].abs(),
                                    ))
                                .abs() /
                            max(SMLNUM, XNORM.value) /
                            EPS;
                    if (INFO.value != 0 && INFO.value != 1) {
                      RES += ONE / EPS;
                    }
                    KNT.value = KNT.value + 1;
                    if (RES > RMAX.value) {
                      LMAX.value = KNT.value;
                      RMAX.value = RES;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
