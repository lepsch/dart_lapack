import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlaexc.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

import 'dhst01.dart';

void dget34(
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
  const TWO = 2.0, THREE = 3.0;
  const LWORK = 32;
  int I,
      IA,
      IA11,
      IA12,
      IA21,
      IA22,
      IAM,
      IB,
      IC,
      IC11,
      IC12,
      IC21,
      IC22,
      ICM,
      J;
  double BIGNUM, EPS, RES, SMLNUM, TNRM;
  final Q = Matrix<double>(4, 4),
      T = Matrix<double>(4, 4),
      T1 = Matrix<double>(4, 4);
  final RESULT = Array<double>(2),
      VAL = Array<double>(9),
      VM = Array<double>(2),
      WORK = Array<double>(LWORK);
  final INFO = Box(0);

  // Get machine parameters

  EPS = dlamch('P');
  SMLNUM = dlamch('S') / EPS;
  BIGNUM = ONE / SMLNUM;

  // Set up test case parameters

  VAL[1] = ZERO;
  VAL[2] = sqrt(SMLNUM);
  VAL[3] = ONE;
  VAL[4] = TWO;
  VAL[5] = sqrt(BIGNUM);
  VAL[6] = -sqrt(SMLNUM);
  VAL[7] = -ONE;
  VAL[8] = -TWO;
  VAL[9] = -sqrt(BIGNUM);
  VM[1] = ONE;
  VM[2] = ONE + TWO * EPS;
  dcopy(16, VAL(4), 0, T(1, 1).asArray(), 1);

  NINFO[1] = 0;
  NINFO[2] = 0;
  KNT.value = 0;
  LMAX.value = 0;
  RMAX.value = ZERO;

  // Begin test loop

  for (IA = 1; IA <= 9; IA++) {
    for (IAM = 1; IAM <= 2; IAM++) {
      for (IB = 1; IB <= 9; IB++) {
        for (IC = 1; IC <= 9; IC++) {
          T[1][1] = VAL[IA] * VM[IAM];
          T[2][2] = VAL[IC];
          T[1][2] = VAL[IB];
          T[2][1] = ZERO;
          TNRM = max(T[1][1].abs(), max(T[2][2].abs(), T[1][2].abs()));
          dcopy(16, T.asArray(), 1, T1.asArray(), 1);
          dcopy(16, VAL(1), 0, Q.asArray(), 1);
          dcopy(4, VAL(3), 0, Q.asArray(), 5);
          dlaexc(true, 2, T, 4, Q, 4, 1, 1, 1, WORK, INFO);
          if (INFO.value != 0) NINFO[INFO.value]++;
          dhst01(2, 1, 2, T1, 4, T, 4, Q, 4, WORK, LWORK, RESULT);
          RES = RESULT[1] + RESULT[2];
          if (INFO.value != 0) RES += ONE / EPS;
          if (T[1][1] != T1[2][2]) RES += ONE / EPS;
          if (T[2][2] != T1[1][1]) RES += ONE / EPS;
          if (T[2][1] != ZERO) RES += ONE / EPS;
          KNT.value++;
          if (RES > RMAX.value) {
            LMAX.value = KNT.value;
            RMAX.value = RES;
          }
        }
      }
    }
  }

  for (IA = 1; IA <= 5; IA++) {
    for (IAM = 1; IAM <= 2; IAM++) {
      for (IB = 1; IB <= 5; IB++) {
        for (IC11 = 1; IC11 <= 5; IC11++) {
          for (IC12 = 2; IC12 <= 5; IC12++) {
            for (IC21 = 2; IC21 <= 4; IC21++) {
              for (IC22 = -1; IC22 <= 1; IC22 += 2) {
                T[1][1] = VAL[IA] * VM[IAM];
                T[1][2] = VAL[IB];
                T[1][3] = -TWO * VAL[IB];
                T[2][1] = ZERO;
                T[2][2] = VAL[IC11];
                T[2][3] = VAL[IC12];
                T[3][1] = ZERO;
                T[3][2] = -VAL[IC21];
                T[3][3] = VAL[IC11] * IC22;
                TNRM = max(
                    T[1][1].abs(),
                    max(
                      T[1][2].abs(),
                      max(
                        T[1][3].abs(),
                        max(
                          T[2][2].abs(),
                          max(
                            T[2][3].abs(),
                            max(T[3][2].abs(), T[3][3].abs()),
                          ),
                        ),
                      ),
                    ));
                dcopy(16, T.asArray(), 1, T1.asArray(), 1);
                dcopy(16, VAL(1), 0, Q.asArray(), 1);
                dcopy(4, VAL(3), 0, Q.asArray(), 5);
                dlaexc(true, 3, T, 4, Q, 4, 1, 1, 2, WORK, INFO);
                if (INFO.value != 0) NINFO[INFO.value]++;
                dhst01(3, 1, 3, T1, 4, T, 4, Q, 4, WORK, LWORK, RESULT);
                RES = RESULT[1] + RESULT[2];
                if (INFO.value == 0) {
                  if (T1[1][1] != T[3][3]) RES += ONE / EPS;
                  if (T[3][1] != ZERO) RES += ONE / EPS;
                  if (T[3][2] != ZERO) RES += ONE / EPS;
                  if (T[2][1] != 0 &&
                      (T[1][1] != T[2][2] ||
                          sign(ONE, T[1][2]) == sign(ONE, T[2][1]))) {
                    RES += ONE / EPS;
                  }
                }
                KNT.value++;
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

  for (IA11 = 1; IA11 <= 5; IA11++) {
    for (IA12 = 2; IA12 <= 5; IA12++) {
      for (IA21 = 2; IA21 <= 4; IA21++) {
        for (IA22 = -1; IA22 <= 1; IA22 += 2) {
          for (ICM = 1; ICM <= 2; ICM++) {
            for (IB = 1; IB <= 5; IB++) {
              for (IC = 1; IC <= 5; IC++) {
                T[1][1] = VAL[IA11];
                T[1][2] = VAL[IA12];
                T[1][3] = -TWO * VAL[IB];
                T[2][1] = -VAL[IA21];
                T[2][2] = VAL[IA11] * IA22;
                T[2][3] = VAL[IB];
                T[3][1] = ZERO;
                T[3][2] = ZERO;
                T[3][3] = VAL[IC] * VM[ICM];
                TNRM = max(
                    T[1][1].abs(),
                    max(
                      T[1][2].abs(),
                      max(
                        T[1][3].abs(),
                        max(
                          T[2][2].abs(),
                          max(
                            T[2][3].abs(),
                            max(T[3][2].abs(), T[3][3].abs()),
                          ),
                        ),
                      ),
                    ));
                dcopy(16, T.asArray(), 1, T1.asArray(), 1);
                dcopy(16, VAL(1), 0, Q.asArray(), 1);
                dcopy(4, VAL(3), 0, Q.asArray(), 5);
                dlaexc(true, 3, T, 4, Q, 4, 1, 2, 1, WORK, INFO);
                if (INFO.value != 0) NINFO[INFO.value]++;
                dhst01(3, 1, 3, T1, 4, T, 4, Q, 4, WORK, LWORK, RESULT);
                RES = RESULT[1] + RESULT[2];
                if (INFO.value == 0) {
                  if (T1[3][3] != T[1][1]) RES += ONE / EPS;
                  if (T[2][1] != ZERO) RES += ONE / EPS;
                  if (T[3][1] != ZERO) RES += ONE / EPS;
                  if (T[3][2] != 0 &&
                      (T[2][2] != T[3][3] ||
                          sign(ONE, T[2][3]) == sign(ONE, T[3][2]))) {
                    RES += ONE / EPS;
                  }
                }
                KNT.value++;
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

  for (IA11 = 1; IA11 <= 5; IA11++) {
    for (IA12 = 2; IA12 <= 5; IA12++) {
      for (IA21 = 2; IA21 <= 4; IA21++) {
        for (IA22 = -1; IA22 <= 1; IA22 += 2) {
          for (IB = 1; IB <= 5; IB++) {
            for (IC11 = 3; IC11 <= 4; IC11++) {
              for (IC12 = 3; IC12 <= 4; IC12++) {
                for (IC21 = 3; IC21 <= 4; IC21++) {
                  for (IC22 = -1; IC22 <= 1; IC22 += 2) {
                    for (ICM = 5; ICM <= 7; ICM++) {
                      IAM = 1;
                      T[1][1] = VAL[IA11] * VM[IAM];
                      T[1][2] = VAL[IA12] * VM[IAM];
                      T[1][3] = -TWO * VAL[IB];
                      T[1][4] = HALF * VAL[IB];
                      T[2][1] = -T[1][2] * VAL[IA21];
                      T[2][2] = VAL[IA11] * IA22 * VM[IAM];
                      T[2][3] = VAL[IB];
                      T[2][4] = THREE * VAL[IB];
                      T[3][1] = ZERO;
                      T[3][2] = ZERO;
                      T[3][3] = VAL[IC11] * VAL[ICM].abs();
                      T[3][4] = VAL[IC12] * VAL[ICM].abs();
                      T[4][1] = ZERO;
                      T[4][2] = ZERO;
                      T[4][3] = -T[3][4] * VAL[IC21] * VAL[ICM].abs();
                      T[4][4] = VAL[IC11] * IC22 * VAL[ICM].abs();
                      TNRM = ZERO;
                      for (I = 1; I <= 4; I++) {
                        for (J = 1; J <= 4; J++) {
                          TNRM = max(TNRM, T[I][J].abs());
                        }
                      }
                      dcopy(16, T.asArray(), 1, T1.asArray(), 1);
                      dcopy(16, VAL(1), 0, Q.asArray(), 1);
                      dcopy(4, VAL(3), 0, Q.asArray(), 5);
                      dlaexc(true, 4, T, 4, Q, 4, 1, 2, 2, WORK, INFO);
                      if (INFO.value != 0) {
                        NINFO[INFO.value]++;
                      }
                      dhst01(4, 1, 4, T1, 4, T, 4, Q, 4, WORK, LWORK, RESULT);
                      RES = RESULT[1] + RESULT[2];
                      if (INFO.value == 0) {
                        if (T[3][1] != ZERO) RES += ONE / EPS;
                        if (T[4][1] != ZERO) RES += ONE / EPS;
                        if (T[3][2] != ZERO) RES += ONE / EPS;
                        if (T[4][2] != ZERO) RES += ONE / EPS;
                        if (T[2][1] != 0 &&
                            (T[1][1] != T[2][2] ||
                                sign(ONE, T[1][2]) == sign(ONE, T[2][1]))) {
                          RES += ONE / EPS;
                        }
                        if (T[4][3] != 0 &&
                            (T[3][3] != T[4][4] ||
                                sign(ONE, T[3][4]) == sign(ONE, T[4][3]))) {
                          RES += ONE / EPS;
                        }
                      }
                      KNT.value++;
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
