import 'dart:math';

import 'package:collection/collection.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlasy2.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dget32(
  final Box<double> RMAX,
  final Box<int> LMAX,
  final Box<int> NINFO,
  final Box<int> KNT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  const TWO = 2.0, FOUR = 4.0, EIGHT = 8.0;
  bool LTRANL, LTRANR;
  int IB, IB1, IB2, IB3, ISGN, ITL, ITLSCL, ITR, ITRANL, ITRANR, ITRSCL, N1, N2;
  double BIGNUM, DEN, EPS, RES, SGN, SMLNUM, TMP, TNRM, XNRM;
  final B = Matrix<double>(2, 2),
      TL = Matrix<double>(2, 2),
      TR = Matrix<double>(2, 2),
      X = Matrix<double>(2, 2);
  final VAL = Array<double>(3);
  final INFO = Box(0);
  final XNORM = Box(0.0), SCALE = Box(0.0);
  final ITVAL = Matrix3d.fromSlice(
      Array.fromList([
        [8, 4, 2, 1, 4, 8, 1, 2],
        [2, 1, 8, 4, 1, 2, 4, 8],
        [9, 4, 2, 1, 4, 9, 1, 2],
        [2, 1, 9, 4, 1, 2, 4, 9],
      ].flattened.toList()),
      [2, 2, 8]);

  // Get machine parameters

  EPS = dlamch('P');
  SMLNUM = dlamch('S') / EPS;
  BIGNUM = ONE / SMLNUM;

  // Set up test case parameters

  VAL[1] = sqrt(SMLNUM);
  VAL[2] = ONE;
  VAL[3] = sqrt(BIGNUM);

  KNT.value = 0;
  NINFO.value = 0;
  LMAX.value = 0;
  RMAX.value = ZERO;

  // Begin test loop

  for (ITRANL = 0; ITRANL <= 1; ITRANL++) {
    for (ITRANR = 0; ITRANR <= 1; ITRANR++) {
      for (ISGN = -1; 2 < 0 ? ISGN >= 1 : ISGN <= 1; ISGN += 2) {
        SGN = ISGN.toDouble();
        LTRANL = ITRANL == 1;
        LTRANR = ITRANR == 1;

        N1 = 1;
        N2 = 1;
        for (ITL = 1; ITL <= 3; ITL++) {
          for (ITR = 1; ITR <= 3; ITR++) {
            for (IB = 1; IB <= 3; IB++) {
              TL[1][1] = VAL[ITL];
              TR[1][1] = VAL[ITR];
              B[1][1] = VAL[IB];
              KNT.value = KNT.value + 1;
              dlasy2(LTRANL, LTRANR, ISGN, N1, N2, TL, 2, TR, 2, B, 2, SCALE, X,
                  2, XNORM, INFO);
              if (INFO.value != 0) NINFO.value = NINFO.value + 1;
              RES = ((TL[1][1] + SGN * TR[1][1]) * X[1][1] -
                      SCALE.value * B[1][1])
                  .abs();
              if (INFO.value == 0) {
                DEN = max(
                    EPS * ((TR[1][1].abs() + TL[1][1]).abs() * X[1][1].abs()),
                    SMLNUM);
              } else {
                DEN = SMLNUM * max(X[1][1].abs(), ONE);
              }
              RES = RES / DEN;
              if (SCALE.value > ONE) RES = RES + ONE / EPS;
              RES = RES +
                  (XNORM.value - X[1][1].abs()).abs() /
                      max(SMLNUM, XNORM.value) /
                      EPS;
              if (INFO.value != 0 && INFO.value != 1) RES = RES + ONE / EPS;
              if (RES > RMAX.value) {
                LMAX.value = KNT.value;
                RMAX.value = RES;
              }
            }
          }
        }

        N1 = 2;
        N2 = 1;
        for (ITL = 1; ITL <= 8; ITL++) {
          for (ITLSCL = 1; ITLSCL <= 3; ITLSCL++) {
            for (ITR = 1; ITR <= 3; ITR++) {
              for (IB1 = 1; IB1 <= 3; IB1++) {
                for (IB2 = 1; IB2 <= 3; IB2++) {
                  B[1][1] = VAL[IB1];
                  B[2][1] = -FOUR * VAL[IB2];
                  TL[1][1] = ITVAL[1][1][ITL] * VAL[ITLSCL];
                  TL[2][1] = ITVAL[2][1][ITL] * VAL[ITLSCL];
                  TL[1][2] = ITVAL[1][2][ITL] * VAL[ITLSCL];
                  TL[2][2] = ITVAL[2][2][ITL] * VAL[ITLSCL];
                  TR[1][1] = VAL[ITR];
                  KNT.value = KNT.value + 1;
                  dlasy2(LTRANL, LTRANR, ISGN, N1, N2, TL, 2, TR, 2, B, 2,
                      SCALE, X, 2, XNORM, INFO);
                  if (INFO.value != 0) NINFO.value = NINFO.value + 1;
                  if (LTRANL) {
                    TMP = TL[1][2];
                    TL[1][2] = TL[2][1];
                    TL[2][1] = TMP;
                  }
                  RES = ((TL[1][1] + SGN * TR[1][1]) * X[1][1] +
                          TL[1][2] * X[2][1] -
                          SCALE.value * B[1][1])
                      .abs();
                  RES = RES +
                      ((TL[2][2] + SGN * TR[1][1]) * X[2][1] +
                              TL[2][1] * X[1][1] -
                              SCALE.value * B[2][1])
                          .abs();
                  TNRM = TR[1][1].abs() +
                      TL[1][1].abs() +
                      TL[1][2].abs() +
                      TL[2][1].abs() +
                      TL[2][2].abs();
                  XNRM = max(X[1][1].abs(), X[2][1]).abs();
                  DEN = max(SMLNUM, max(SMLNUM * XNRM, (TNRM * EPS) * XNRM));
                  RES = RES / DEN;
                  if (SCALE.value > ONE) RES = RES + ONE / EPS;
                  RES = RES +
                      (XNORM.value - XNRM).abs() /
                          max(SMLNUM, XNORM.value) /
                          EPS;
                  if (RES > RMAX.value) {
                    LMAX.value = KNT.value;
                    RMAX.value = RES;
                  }
                }
              }
            }
          }
        }

        N1 = 1;
        N2 = 2;
        for (ITR = 1; ITR <= 8; ITR++) {
          for (ITRSCL = 1; ITRSCL <= 3; ITRSCL++) {
            for (ITL = 1; ITL <= 3; ITL++) {
              for (IB1 = 1; IB1 <= 3; IB1++) {
                for (IB2 = 1; IB2 <= 3; IB2++) {
                  B[1][1] = VAL[IB1];
                  B[1][2] = -TWO * VAL[IB2];
                  TR[1][1] = ITVAL[1][1][ITR] * VAL[ITRSCL];
                  TR[2][1] = ITVAL[2][1][ITR] * VAL[ITRSCL];
                  TR[1][2] = ITVAL[1][2][ITR] * VAL[ITRSCL];
                  TR[2][2] = ITVAL[2][2][ITR] * VAL[ITRSCL];
                  TL[1][1] = VAL[ITL];
                  KNT.value = KNT.value + 1;
                  dlasy2(LTRANL, LTRANR, ISGN, N1, N2, TL, 2, TR, 2, B, 2,
                      SCALE, X, 2, XNORM, INFO);
                  if (INFO.value != 0) NINFO.value = NINFO.value + 1;
                  if (LTRANR) {
                    TMP = TR[1][2];
                    TR[1][2] = TR[2][1];
                    TR[2][1] = TMP;
                  }
                  TNRM = TL[1][1].abs() +
                      TR[1][1].abs() +
                      TR[1][2].abs() +
                      TR[2][2].abs() +
                      TR[2][1].abs();
                  XNRM = X[1][1].abs() + X[1][2].abs();
                  RES = (((TL[1][1] + SGN * TR[1][1])) * X[1][1] +
                          (SGN * TR[2][1]) * X[1][2] -
                          (SCALE.value * B[1][1]))
                      .abs();
                  RES = RES +
                      (((TL[1][1] + SGN * TR[2][2])) * X[1][2] +
                              (SGN * TR[1][2]) * X[1][1] -
                              (SCALE.value * B[1][2]))
                          .abs();
                  DEN = max(SMLNUM, max(SMLNUM * XNRM, (TNRM * EPS) * XNRM));
                  RES = RES / DEN;
                  if (SCALE.value > ONE) RES = RES + ONE / EPS;
                  RES = RES +
                      (XNORM.value - XNRM).abs() /
                          max(SMLNUM, XNORM.value) /
                          EPS;
                  if (RES > RMAX.value) {
                    LMAX.value = KNT.value;
                    RMAX.value = RES;
                  }
                }
              }
            }
          }
        }

        N1 = 2;
        N2 = 2;
        for (ITR = 1; ITR <= 8; ITR++) {
          for (ITRSCL = 1; ITRSCL <= 3; ITRSCL++) {
            for (ITL = 1; ITL <= 8; ITL++) {
              for (ITLSCL = 1; ITLSCL <= 3; ITLSCL++) {
                for (IB1 = 1; IB1 <= 3; IB1++) {
                  for (IB2 = 1; IB2 <= 3; IB2++) {
                    for (IB3 = 1; IB3 <= 3; IB3++) {
                      B[1][1] = VAL[IB1];
                      B[2][1] = -FOUR * VAL[IB2];
                      B[1][2] = -TWO * VAL[IB3];
                      B[2][2] = EIGHT * min(VAL[IB1], min(VAL[IB2], VAL[IB3]));
                      TR[1][1] = ITVAL[1][1][ITR] * VAL[ITRSCL];
                      TR[2][1] = ITVAL[2][1][ITR] * VAL[ITRSCL];
                      TR[1][2] = ITVAL[1][2][ITR] * VAL[ITRSCL];
                      TR[2][2] = ITVAL[2][2][ITR] * VAL[ITRSCL];
                      TL[1][1] = ITVAL[1][1][ITL] * VAL[ITLSCL];
                      TL[2][1] = ITVAL[2][1][ITL] * VAL[ITLSCL];
                      TL[1][2] = ITVAL[1][2][ITL] * VAL[ITLSCL];
                      TL[2][2] = ITVAL[2][2][ITL] * VAL[ITLSCL];
                      KNT.value = KNT.value + 1;
                      dlasy2(LTRANL, LTRANR, ISGN, N1, N2, TL, 2, TR, 2, B, 2,
                          SCALE, X, 2, XNORM, INFO);
                      if (INFO.value != 0) NINFO.value = NINFO.value + 1;
                      if (LTRANR) {
                        TMP = TR[1][2];
                        TR[1][2] = TR[2][1];
                        TR[2][1] = TMP;
                      }
                      if (LTRANL) {
                        TMP = TL[1][2];
                        TL[1][2] = TL[2][1];
                        TL[2][1] = TMP;
                      }
                      TNRM = TR[1][1].abs() +
                          TR[2][1].abs() +
                          TR[1][2].abs() +
                          TR[2][2].abs() +
                          TL[1][1].abs() +
                          TL[2][1].abs() +
                          TL[1][2].abs() +
                          TL[2][2].abs();
                      XNRM = max(
                        X[1][1].abs() + X[1][2].abs(),
                        X[2][1].abs() + X[2][2],
                      ).abs();
                      RES = (((TL[1][1] + SGN * TR[1][1])) * X[1][1] +
                              (SGN * TR[2][1]) * X[1][2] +
                              TL[1][2] * X[2][1] -
                              (SCALE.value * B[1][1]))
                          .abs();
                      RES = RES +
                          (TL[1][1] * X[1][2] +
                                  (SGN * TR[1][2]) * X[1][1] +
                                  (SGN * TR[2][2]) * X[1][2] +
                                  TL[1][2] * X[2][2] -
                                  (SCALE.value * B[1][2]))
                              .abs();
                      RES = RES +
                          (TL[2][1] * X[1][1] +
                                  (SGN * TR[1][1]) * X[2][1] +
                                  (SGN * TR[2][1]) * X[2][2] +
                                  TL[2][2] * X[2][1] -
                                  (SCALE.value * B[2][1]))
                              .abs();
                      RES = RES +
                          (((TL[2][2] + SGN * TR[2][2])) * X[2][2] +
                                  (SGN * TR[1][2]) * X[2][1] +
                                  TL[2][1] * X[1][2] -
                                  (SCALE.value * B[2][2]))
                              .abs();
                      DEN =
                          max(SMLNUM, max(SMLNUM * XNRM, (TNRM * EPS) * XNRM));
                      RES = RES / DEN;
                      if (SCALE.value > ONE) RES = RES + ONE / EPS;
                      RES = RES +
                          (XNORM.value - XNRM).abs() /
                              max(SMLNUM, XNORM.value) /
                              EPS;
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
}
