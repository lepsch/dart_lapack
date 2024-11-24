// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlapy2.dart';
import 'package:dart_lapack/src/intrinsics/sign.dart';
import 'package:dart_lapack/src/install/dlamch.dart';

void dlanv2(
  final Box<double> A,
  final Box<double> B,
  final Box<double> C,
  final Box<double> D,
  final Box<double> RT1R,
  final Box<double> RT1I,
  final Box<double> RT2R,
  final Box<double> RT2I,
  final Box<double> CS,
  final Box<double> SN,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, HALF = 0.5, ONE = 1.0, TWO = 2.0;
  const MULTPL = 4.0;
  double AA,
      BB,
      BCMAX,
      BCMIS,
      CC,
      CS1,
      DD,
      EPS,
      P,
      SAB,
      SAC,
      SCALE,
      SIGMA = 0,
      SN1,
      TAU,
      TEMP,
      Z,
      SAFMIN,
      SAFMN2,
      SAFMX2;
  int COUNT = 0;

  SAFMIN = dlamch('S');
  EPS = dlamch('P');
  SAFMN2 = pow(
    dlamch('B'),
    (log(SAFMIN / EPS) / log(dlamch('B')) ~/ TWO),
  ).toDouble();
  SAFMX2 = ONE / SAFMN2;
  if (C.value == ZERO) {
    CS.value = ONE;
    SN.value = ZERO;
  } else if (B.value == ZERO) {
    // Swap rows and columns

    CS.value = ZERO;
    SN.value = ONE;
    TEMP = D.value;
    D.value = A.value;
    A.value = TEMP;
    B.value = -C.value;
    C.value = ZERO;
  } else if ((A.value - D.value) == ZERO &&
      sign(ONE, B.value) != sign(ONE, C.value)) {
    CS.value = ONE;
    SN.value = ZERO;
  } else {
    TEMP = A.value - D.value;
    P = HALF * TEMP;
    BCMAX = max(B.value.abs(), C.value.abs());
    BCMIS = min(B.value.abs(), C.value.abs()) *
        sign(ONE, B.value) *
        sign(ONE, C.value);
    SCALE = max(P.abs(), BCMAX);
    Z = (P / SCALE) * P + (BCMAX / SCALE) * BCMIS;

    // If Z is of the order of the machine accuracy, postpone the
    // decision on the nature of eigenvalues

    if (Z >= MULTPL * EPS) {
      // Real eigenvalues. Compute A and D.

      Z = P + sign(sqrt(SCALE) * sqrt(Z), P);
      A.value = D.value + Z;
      D.value -= (BCMAX / Z) * BCMIS;

      // Compute B and the rotation matrix

      TAU = dlapy2(C.value, Z);
      CS.value = Z / TAU;
      SN.value = C.value / TAU;
      B.value -= C.value;
      C.value = ZERO;
    } else {
      // Complex eigenvalues, or real (almost) equal eigenvalues.
      // Make diagonal elements equal.

      COUNT = 0;
      SIGMA = B.value + C.value;
      while (true) {
        COUNT++;
        SCALE = max(TEMP.abs(), SIGMA.abs());
        if (SCALE >= SAFMX2) {
          SIGMA *= SAFMN2;
          TEMP *= SAFMN2;
          if (COUNT <= 20) continue;
        }
        if (SCALE <= SAFMN2) {
          SIGMA *= SAFMX2;
          TEMP *= SAFMX2;
          if (COUNT <= 20) continue;
        }
        break;
      }
      P = HALF * TEMP;
      TAU = dlapy2(SIGMA, TEMP);
      CS.value = sqrt(HALF * (ONE + SIGMA.abs() / TAU));
      SN.value = -(P / (TAU * CS.value)) * sign(ONE, SIGMA);

      // Compute [ AA  BB ] = [ A  B ] [ CS -SN ]
      // [ CC  DD ]   [ C  D ] [ SN  CS ]

      AA = A.value * CS.value + B.value * SN.value;
      BB = -A.value * SN.value + B.value * CS.value;
      CC = C.value * CS.value + D.value * SN.value;
      DD = -C.value * SN.value + D.value * CS.value;

      // Compute [ A  B ] = [ CS  SN ] [ AA  BB ]
      // [ C  D ]   [-SN  CS ] [ CC  DD ]

      A.value = AA * CS.value + CC * SN.value;
      B.value = BB * CS.value + DD * SN.value;
      C.value = -AA * SN.value + CC * CS.value;
      D.value = -BB * SN.value + DD * CS.value;

      TEMP = HALF * (A.value + D.value);
      A.value = TEMP;
      D.value = TEMP;

      if (C.value != ZERO) {
        if (B.value != ZERO) {
          if (sign(ONE, B.value) == sign(ONE, C.value)) {
            // Real eigenvalues: reduce to upper triangular form

            SAB = sqrt(B.value.abs());
            SAC = sqrt(C.value.abs());
            P = sign(SAB * SAC, C.value);
            TAU = ONE / sqrt((B.value + C.value).abs());
            A.value = TEMP + P;
            D.value = TEMP - P;
            B.value -= C.value;
            C.value = ZERO;
            CS1 = SAB * TAU;
            SN1 = SAC * TAU;
            TEMP = CS.value * CS1 - SN.value * SN1;
            SN.value = CS.value * SN1 + SN.value * CS1;
            CS.value = TEMP;
          }
        } else {
          B.value = -C.value;
          C.value = ZERO;
          TEMP = CS.value;
          CS.value = -SN.value;
          SN.value = TEMP;
        }
      }
    }
  }

  // Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).

  RT1R.value = A.value;
  RT2R.value = D.value;
  if (C.value == ZERO) {
    RT1I.value = ZERO;
    RT2I.value = ZERO;
  } else {
    RT1I.value = sqrt(B.value.abs()) * sqrt(C.value.abs());
    RT2I.value = -RT1I.value;
  }
}
