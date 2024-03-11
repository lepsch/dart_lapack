import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/install/dlamch.dart';

void dlartgp(
  final double F,
  final double G,
  final Box<double> CS,
  final Box<double> SN,
  final Box<double> R,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;
  const ONE = 1.0;
  const TWO = 2.0;
  // bool              FIRST
  int COUNT = 0, I;
  double EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE;

  // IF( FIRST ) THEN
  SAFMIN = dlamch('S');
  EPS = dlamch('E');
  SAFMN2 = pow(dlamch('B'), (log(SAFMIN / EPS) / log(dlamch('B')) ~/ TWO))
      .toDouble();
  SAFMX2 = ONE / SAFMN2;
  // FIRST = false;
  // END IF
  if (G == ZERO) {
    CS.value = sign(ONE, F).toDouble();
    SN.value = ZERO;
    R.value = (F).abs();
  } else if (F == ZERO) {
    CS.value = ZERO;
    SN.value = sign(ONE, G).toDouble();
    R.value = (G).abs();
  } else {
    F1 = F;
    G1 = G;
    SCALE = max((F1).abs(), (G1).abs());
    if (SCALE >= SAFMX2) {
      COUNT = 0;
      do {
        COUNT += 1;
        F1 = F1 * SAFMN2;
        G1 = G1 * SAFMN2;
        SCALE = max((F1).abs(), (G1).abs());
      } while (SCALE >= SAFMX2 && COUNT < 20);
      R.value = sqrt(pow(F1, 2) + pow(G1, 2));
      CS.value = F1 / R.value;
      SN.value = G1 / R.value;
      for (I = 1; I <= COUNT; I++) {
        R.value = R.value * SAFMX2;
      }
    } else if (SCALE <= SAFMN2) {
      COUNT = 0;
      do {
        COUNT++;
        F1 = F1 * SAFMX2;
        G1 = G1 * SAFMX2;
        SCALE = max((F1).abs(), (G1).abs());
      } while (SCALE <= SAFMN2);
      R.value = sqrt(pow(F1, 2) + pow(G1, 2));
      CS.value = F1 / R.value;
      SN.value = G1 / R.value;
      for (I = 1; I <= COUNT; I++) {
        R.value = R.value * SAFMN2;
      }
    } else {
      R.value = sqrt(pow(F1, 2) + pow(G1, 2));
      CS.value = F1 / R.value;
      SN.value = G1 / R.value;
    }
    if (R.value < ZERO) {
      CS.value = -CS.value;
      SN.value = -SN.value;
      R.value = -R.value;
    }
  }
}
