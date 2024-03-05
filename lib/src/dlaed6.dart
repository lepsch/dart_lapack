import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dlaed6(
  final int KNITER,
  final bool ORGATI,
  final double RHO,
  final Array<double> D_,
  final Array<double> Z_,
  final double FINIT,
  final Box<double> TAU,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final Z = Z_.having();
  const MAXIT = 40;
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, THREE = 3.0, FOUR = 4.0, EIGHT = 8.0;
  final DSCALE = Array<double>(3), ZSCALE = Array<double>(3);
  bool SCALE;
  int I, ITER, NITER;
  double A,
      B,
      BASE,
      C,
      DDF,
      DF,
      EPS,
      ERRETM,
      ETA,
      F,
      FC,
      SCLFAC,
      SCLINV = 0,
      SMALL1,
      SMALL2,
      SMINV1,
      SMINV2,
      TEMP,
      TEMP1,
      TEMP2,
      TEMP3,
      TEMP4,
      LBD,
      UBD;

  INFO.value = 0;

  if (ORGATI) {
    LBD = D[2];
    UBD = D[3];
  } else {
    LBD = D[1];
    UBD = D[2];
  }
  if (FINIT < ZERO) {
    LBD = ZERO;
  } else {
    UBD = ZERO;
  }

  NITER = 1;
  TAU.value = ZERO;
  if (KNITER == 2) {
    if (ORGATI) {
      TEMP = (D[3] - D[2]) / TWO;
      C = RHO + Z[1] / ((D[1] - D[2]) - TEMP);
      A = C * (D[2] + D[3]) + Z[2] + Z[3];
      B = C * D[2] * D[3] + Z[2] * D[3] + Z[3] * D[2];
    } else {
      TEMP = (D[1] - D[2]) / TWO;
      C = RHO + Z[3] / ((D[3] - D[2]) - TEMP);
      A = C * (D[1] + D[2]) + Z[1] + Z[2];
      B = C * D[1] * D[2] + Z[1] * D[2] + Z[2] * D[1];
    }
    TEMP = max(A.abs(), max(B.abs(), C.abs()));
    A = A / TEMP;
    B = B / TEMP;
    C = C / TEMP;
    if (C == ZERO) {
      TAU.value = B / A;
    } else if (A <= ZERO) {
      TAU.value = (A - sqrt((A * A - FOUR * B * C).abs())) / (TWO * C);
    } else {
      TAU.value = TWO * B / (A + sqrt((A * A - FOUR * B * C).abs()));
    }
    if (TAU.value < LBD || TAU.value > UBD) TAU.value = (LBD + UBD) / TWO;
    if (D[1] == TAU.value || D[2] == TAU.value || D[3] == TAU.value) {
      TAU.value = ZERO;
    } else {
      TEMP = FINIT +
          TAU.value * Z[1] / (D[1] * (D[1] - TAU.value)) +
          TAU.value * Z[2] / (D[2] * (D[2] - TAU.value)) +
          TAU.value * Z[3] / (D[3] * (D[3] - TAU.value));
      if (TEMP <= ZERO) {
        LBD = TAU.value;
      } else {
        UBD = TAU.value;
      }
      if (FINIT.abs() <= TEMP.abs()) TAU.value = ZERO;
    }
  }

  // get machine parameters for possible scaling to avoid overflow

  // modified by Sven: parameters SMALL1, SMINV1, SMALL2,
  // SMINV2, EPS are not SAVEd anymore between one call to the
  // others but recomputed at each call

  EPS = dlamch('Epsilon');
  BASE = dlamch('Base');
  SMALL1 = pow(BASE, log(dlamch('SafMin')) / log(BASE) ~/ THREE).toDouble();
  SMINV1 = ONE / SMALL1;
  SMALL2 = SMALL1 * SMALL1;
  SMINV2 = SMINV1 * SMINV1;

  // Determine if scaling of inputs necessary to avoid overflow
  // when computing 1/TEMP**3

  if (ORGATI) {
    TEMP = min((D[2] - TAU.value).abs(), (D[3] - TAU.value).abs());
  } else {
    TEMP = min((D[1] - TAU.value).abs(), (D[2] - TAU.value).abs());
  }
  SCALE = false;
  if (TEMP <= SMALL1) {
    SCALE = true;
    if (TEMP <= SMALL2) {
      // Scale up by power of radix nearest 1/SAFMIN**(2/3)

      SCLFAC = SMINV2;
      SCLINV = SMALL2;
    } else {
      // Scale up by power of radix nearest 1/SAFMIN**(1/3)

      SCLFAC = SMINV1;
      SCLINV = SMALL1;
    }

    // Scaling up safe because D, Z, TAU.value scaled elsewhere to be O(1)

    for (I = 1; I <= 3; I++) {
      DSCALE[I] = D[I] * SCLFAC;
      ZSCALE[I] = Z[I] * SCLFAC;
    }
    TAU.value = TAU.value * SCLFAC;
    LBD = LBD * SCLFAC;
    UBD = UBD * SCLFAC;
  } else {
    // Copy D and Z to DSCALE and ZSCALE

    for (I = 1; I <= 3; I++) {
      DSCALE[I] = D[I];
      ZSCALE[I] = Z[I];
    }
  }

  FC = ZERO;
  DF = ZERO;
  DDF = ZERO;
  for (I = 1; I <= 3; I++) {
    TEMP = ONE / (DSCALE[I] - TAU.value);
    TEMP1 = ZSCALE[I] * TEMP;
    TEMP2 = TEMP1 * TEMP;
    TEMP3 = TEMP2 * TEMP;
    FC = FC + TEMP1 / DSCALE[I];
    DF = DF + TEMP2;
    DDF = DDF + TEMP3;
  }
  F = FINIT + TAU.value * FC;

  done:
  while (true) {
    if ((F).abs() <= ZERO) break done;
    if (F <= ZERO) {
      LBD = TAU.value;
    } else {
      UBD = TAU.value;
    }

    // Iteration begins -- Use Gragg-Thornton-Warner cubic convergent
    // scheme

    // It is not hard to see that

    // 1) Iterations will go up monotonically
    // if FINIT < 0;

    // 2) Iterations will go down monotonically
    // if FINIT > 0.

    ITER = NITER + 1;

    for (NITER = ITER; NITER <= MAXIT; NITER++) {
      if (ORGATI) {
        TEMP1 = DSCALE[2] - TAU.value;
        TEMP2 = DSCALE[3] - TAU.value;
      } else {
        TEMP1 = DSCALE[1] - TAU.value;
        TEMP2 = DSCALE[2] - TAU.value;
      }
      A = (TEMP1 + TEMP2) * F - TEMP1 * TEMP2 * DF;
      B = TEMP1 * TEMP2 * F;
      C = F - (TEMP1 + TEMP2) * DF + TEMP1 * TEMP2 * DDF;
      TEMP = max(A.abs(), max(B.abs(), C.abs()));
      A = A / TEMP;
      B = B / TEMP;
      C = C / TEMP;
      if (C == ZERO) {
        ETA = B / A;
      } else if (A <= ZERO) {
        ETA = (A - sqrt((A * A - FOUR * B * C).abs())) / (TWO * C);
      } else {
        ETA = TWO * B / (A + sqrt((A * A - FOUR * B * C).abs()));
      }
      if (F * ETA >= ZERO) {
        ETA = -F / DF;
      }

      TAU.value = TAU.value + ETA;
      if (TAU.value < LBD || TAU.value > UBD) TAU.value = (LBD + UBD) / TWO;

      FC = ZERO;
      ERRETM = ZERO;
      DF = ZERO;
      DDF = ZERO;
      for (I = 1; I <= 3; I++) {
        if ((DSCALE[I] - TAU.value) != ZERO) {
          TEMP = ONE / (DSCALE[I] - TAU.value);
          TEMP1 = ZSCALE[I] * TEMP;
          TEMP2 = TEMP1 * TEMP;
          TEMP3 = TEMP2 * TEMP;
          TEMP4 = TEMP1 / DSCALE[I];
          FC = FC + TEMP4;
          ERRETM = ERRETM + (TEMP4).abs();
          DF = DF + TEMP2;
          DDF = DDF + TEMP3;
        } else {
          break done;
        }
      }
      F = FINIT + TAU.value * FC;
      ERRETM = EIGHT * ((FINIT).abs() + (TAU.value).abs() * ERRETM) +
          (TAU.value).abs() * DF;
      if (((F).abs() <= FOUR * EPS * ERRETM) ||
          ((UBD - LBD) <= FOUR * EPS * (TAU.value).abs())) break done;
      if (F <= ZERO) {
        LBD = TAU.value;
      } else {
        UBD = TAU.value;
      }
    }
    INFO.value = 1;
    break;
  }

  // Undo scaling
  if (SCALE) TAU.value = TAU.value * SCLINV;
}
