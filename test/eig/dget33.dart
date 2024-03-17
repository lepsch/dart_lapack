import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlanv2.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dget33(
  final Box<double> RMAX,
  final Box<int> LMAX,
  final Box<int> NINFO,
  final Box<int> KNT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  const TWO = 2.0, FOUR = 4.0;
  int I1, I2, I3, I4, IM1, IM2, IM3, IM4, J1, J2, J3;
  double BIGNUM, EPS, RES, SMLNUM, SUM, TNRM;
  final Q = Matrix<double>(2, 2),
      T = Matrix<double>(2, 2),
      T1 = Matrix<double>(2, 2),
      T2 = Matrix<double>(2, 2);

  final VAL = Array<double>(4), VM = Array<double>(3);
  final CS = Box(0.0),
      SN = Box(0.0),
      WI1 = Box(0.0),
      WI2 = Box(0.0),
      WR1 = Box(0.0),
      WR2 = Box(0.0);

  // Get machine parameters

  EPS = dlamch('P');
  SMLNUM = dlamch('S') / EPS;
  BIGNUM = ONE / SMLNUM;

  // Set up test case parameters

  VAL[1] = ONE;
  VAL[2] = ONE + TWO * EPS;
  VAL[3] = TWO;
  VAL[4] = TWO - FOUR * EPS;
  VM[1] = SMLNUM;
  VM[2] = ONE;
  VM[3] = BIGNUM;

  KNT.value = 0;
  NINFO.value = 0;
  LMAX.value = 0;
  RMAX.value = ZERO;

  // Begin test loop

  for (I1 = 1; I1 <= 4; I1++) {
    for (I2 = 1; I2 <= 4; I2++) {
      for (I3 = 1; I3 <= 4; I3++) {
        for (I4 = 1; I4 <= 4; I4++) {
          for (IM1 = 1; IM1 <= 3; IM1++) {
            for (IM2 = 1; IM2 <= 3; IM2++) {
              for (IM3 = 1; IM3 <= 3; IM3++) {
                for (IM4 = 1; IM4 <= 3; IM4++) {
                  T[1][1] = VAL[I1] * VM[IM1];
                  T[1][2] = VAL[I2] * VM[IM2];
                  T[2][1] = -VAL[I3] * VM[IM3];
                  T[2][2] = VAL[I4] * VM[IM4];
                  TNRM = max(
                      (T[1][1]).abs(),
                      max(
                        (T[1][2]).abs(),
                        max((T[2][1]).abs(), (T[2][2]).abs()),
                      ));
                  T1[1][1] = T[1][1];
                  T1[1][2] = T[1][2];
                  T1[2][1] = T[2][1];
                  T1[2][2] = T[2][2];
                  Q[1][1] = ONE;
                  Q[1][2] = ZERO;
                  Q[2][1] = ZERO;
                  Q[2][2] = ONE;

                  dlanv2(T.box(1, 1), T.box(1, 2), T.box(2, 1), T.box(2, 2),
                      WR1, WI1, WR2, WI2, CS, SN);
                  for (J1 = 1; J1 <= 2; J1++) {
                    RES = Q[J1][1] * CS.value + Q[J1][2] * SN.value;
                    Q[J1][2] = -Q[J1][1] * SN.value + Q[J1][2] * CS.value;
                    Q[J1][1] = RES;
                  }

                  RES = ZERO;
                  RES = RES +
                      (pow(Q[1][1], 2) + pow(Q[1][2], 2) - ONE).abs() / EPS;
                  RES = RES +
                      (pow(Q[2][2], 2) + pow(Q[2][1], 2) - ONE).abs() / EPS;
                  RES =
                      RES + (Q[1][1] * Q[2][1] + Q[1][2] * Q[2][2]).abs() / EPS;
                  for (J1 = 1; J1 <= 2; J1++) {
                    for (J2 = 1; J2 <= 2; J2++) {
                      T2[J1][J2] = ZERO;
                      for (J3 = 1; J3 <= 2; J3++) {
                        T2[J1][J2] += T1[J1][J3] * Q[J3][J2];
                      }
                    }
                  }
                  for (J1 = 1; J1 <= 2; J1++) {
                    for (J2 = 1; J2 <= 2; J2++) {
                      SUM = T[J1][J2];
                      for (J3 = 1; J3 <= 2; J3++) {
                        SUM -= Q[J3][J1] * T2[J3][J2];
                      }
                      RES += (SUM).abs() / EPS / TNRM;
                    }
                  }
                  if (T[2][1] != ZERO &&
                      (T[1][1] != T[2][2] ||
                          sign(ONE, T[1][2]) * sign(ONE, T[2][1]) > ZERO)) {
                    RES += ONE / EPS;
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
