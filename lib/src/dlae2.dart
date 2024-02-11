import 'dart:math';

import 'package:lapack/src/box.dart';

void dlae2(
  final double A,
  final double B,
  final double C,
  final Box<double> RT1,
  final Box<double> RT2,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0;
  const TWO = 2.0;
  const ZERO = 0.0;
  const HALF = 0.5;
  double AB, ACMN, ACMX, ADF, DF, RT, SM, TB;

  // Compute the eigenvalues

  SM = A + C;
  DF = A - C;
  ADF = (DF).abs();
  TB = B + B;
  AB = (TB).abs();
  if ((A).abs() > (C).abs()) {
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

    // Order of execution important.
    // To get fully accurate smaller eigenvalue,
    // next line needs to be executed in higher precision.

    RT2.value = (ACMX / RT1.value) * ACMN - (B / RT1.value) * B;
  } else if (SM > ZERO) {
    RT1.value = HALF * (SM + RT);

    // Order of execution important.
    // To get fully accurate smaller eigenvalue,
    // next line needs to be executed in higher precision.

    RT2.value = (ACMX / RT1.value) * ACMN - (B / RT1.value) * B;
  } else {
    // Includes case RT1.value = RT2.value = 0

    RT1.value = HALF * RT;
    RT2.value = -HALF * RT;
  }
}
