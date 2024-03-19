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
  const BS = 2.0, HALF = 0.5, TWO = 2.0;

  final AA = Box(A);
  final BB = Box(B);
  var CC = C;
  var DD = D;
  var AB = max(A.abs(), B.abs());
  var CD = max(C.abs(), D.abs());
  var S = 1.0;

  final OV = dlamch('Overflow threshold');
  final UN = dlamch('Safe minimum');
  final EPS = dlamch('Epsilon');
  final BE = BS / (EPS * EPS);

  if (AB >= HALF * OV) {
    AA.value *= HALF;
    BB.value *= HALF;
    S *= TWO;
  }
  if (CD >= HALF * OV) {
    CC *= HALF;
    DD *= HALF;
    S *= HALF;
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
  if (D.abs() <= C.abs()) {
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

  final R = D / C;
  final T = ONE / (C + D * R);
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

  if (R != ZERO) {
    final BR = B * R;
    if (BR != ZERO) {
      return (A + BR) * T;
    }

    return A * T + (B * T) * R;
  }

  return (A + D * (B / C)) * T;
}
