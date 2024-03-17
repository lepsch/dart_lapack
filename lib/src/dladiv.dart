import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/dlamch.dart';

void dladiv(
  final double A,
  final double B,
  final double C,
  final double D,
  final Box<double> P,
  final Box<double> Q,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const BS = 2.0;
  const HALF = 0.5;
  const TWO = 2.0;
  double CC, DD, AB, CD, S, OV, UN, BE, EPS;
  final AA = Box(0.0), BB = Box(0.0);

  AA.value = A;
  BB.value = B;
  CC = C;
  DD = D;
  AB = max((A).abs(), (B).abs());
  CD = max((C).abs(), (D).abs());
  S = 1.0;

  OV = dlamch('Overflow threshold');
  UN = dlamch('Safe minimum');
  EPS = dlamch('Epsilon');
  BE = BS / (EPS * EPS);

  if (AB >= HALF * OV) {
    AA.value = HALF * AA.value;
    BB.value = HALF * BB.value;
    S = TWO * S;
  }
  if (CD >= HALF * OV) {
    CC = HALF * CC;
    DD = HALF * DD;
    S = HALF * S;
  }
  if (AB <= UN * BS / EPS) {
    AA.value *= BE;
    BB.value *= BE;
    S /= BE;
  }
  if (CD <= UN * BS / EPS) {
    CC *= BE;
    DD *= BE;
    S *= BE;
  }
  if ((D).abs() <= (C).abs()) {
    dladiv1(AA, BB.value, CC, DD, P, Q);
  } else {
    dladiv1(BB, AA.value, DD, CC, P, Q);
    Q.value = -Q.value;
  }
  P.value *= S;
  Q.value *= S;
}

void dladiv1(
  final Box<double> A,
  final double B,
  final double C,
  final double D,
  final Box<double> P,
  final Box<double> Q,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0;
  double R, T;

  R = D / C;
  T = ONE / (C + D * R);
  P.value = dladiv2(A.value, B, C, D, R, T);
  A.value = -A.value;
  Q.value = dladiv2(B, A.value, C, D, R, T);
}

double dladiv2(
  final double A,
  final double B,
  final double C,
  final double D,
  final double R,
  final double T,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;
  double BR;

  if (R != ZERO) {
    BR = B * R;
    if (BR != ZERO) {
      return (A + BR) * T;
    } else {
      return A * T + (B * T) * R;
    }
  } else {
    return (A + D * (B / C)) * T;
  }
}
