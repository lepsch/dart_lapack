import 'dart:math';

import 'package:lapack/src/box.dart';

void dlaev2(
  final double A,
  final double B,
  final double C,
  final Box<double> RT1,
  final Box<double> RT2,
  final Box<double> CS1,
  final Box<double> SN1,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0;
  const TWO = 2.0;
  const ZERO = 0.0;
  const HALF = 0.5;
  int SGN1, SGN2;
  double AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM, TB, TN;

  // Compute the eigenvalues

  SM = A + C;
  DF = A - C;
  ADF = DF.abs();
  TB = B + B;
  AB = TB.abs();
  if (A.abs() > C.abs()) {
    ACMX = A;
    ACMN = C;
  } else {
    ACMX = C;
    ACMN = A;
  }
  if (ADF > AB) {
    RT = ADF * sqrt(ONE + pow(AB / ADF, 2));
  } else if (ADF < AB) {
    RT = AB * sqrt(ONE + pow(ADF / AB, 2));
  } else {
    // Includes case AB=ADF=0

    RT = AB * sqrt(TWO);
  }
  if (SM < ZERO) {
    RT1.value = HALF * (SM - RT);
    SGN1 = -1;

    // Order of execution important.
    // To get fully accurate smaller eigenvalue,
    // next line needs to be executed in higher precision.

    RT2.value = (ACMX / RT1.value) * ACMN - (B / RT1.value) * B;
  } else if (SM > ZERO) {
    RT1.value = HALF * (SM + RT);
    SGN1 = 1;

    // Order of execution important.
    // To get fully accurate smaller eigenvalue,
    // next line needs to be executed in higher precision.

    RT2.value = (ACMX / RT1.value) * ACMN - (B / RT1.value) * B;
  } else {
    // Includes case RT1 = RT2 = 0

    RT1.value = HALF * RT;
    RT2.value = -HALF * RT;
    SGN1 = 1;
  }

  // Compute the eigenvector

  if (DF >= ZERO) {
    CS = DF + RT;
    SGN2 = 1;
  } else {
    CS = DF - RT;
    SGN2 = -1;
  }
  ACS = CS.abs();
  if (ACS > AB) {
    CT = -TB / CS;
    SN1.value = ONE / sqrt(ONE + CT * CT);
    CS1.value = CT * SN1.value;
  } else {
    if (AB == ZERO) {
      CS1.value = ONE;
      SN1.value = ZERO;
    } else {
      TN = -CS / TB;
      CS1.value = ONE / sqrt(ONE + TN * TN);
      SN1.value = TN * CS1.value;
    }
  }
  if (SGN1 == SGN2) {
    TN = CS1.value;
    CS1.value = -SN1.value;
    SN1.value = TN;
  }
}
