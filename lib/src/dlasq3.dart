import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/dlasq4.dart';
import 'package:lapack/src/dlasq5.dart';
import 'package:lapack/src/dlasq6.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dlasq3(
  final int I0,
  final Box<int> N0,
  final Array<double> Z_,
  final Box<int> PP,
  final Box<double> DMIN,
  final Box<double> SIGMA,
  final Box<double> DESIG,
  final Box<double> QMAX,
  final Box<int> NFAIL,
  final Box<int> ITER,
  final Box<int> NDIV,
  final bool IEEE,
  final Box<int> TTYPE,
  final Box<double> DMIN1,
  final Box<double> DMIN2,
  final Box<double> DN,
  final Box<double> DN1,
  final Box<double> DN2,
  final Box<double> G,
  final Box<double> TAU,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final Z = Z_.having();
  const CBIAS = 1.50;
  const ZERO = 0.0,
      QURTR = 0.250,
      HALF = 0.5,
      ONE = 1.0,
      TWO = 2.0,
      HUNDRD = 100.0;
  int IPN4, J4, N0IN, NN;
  double EPS, S, T, TEMP, TOL, TOL2;

  N0IN = N0.value;
  EPS = dlamch('Precision');
  TOL = EPS * HUNDRD;
  TOL2 = pow(TOL, 2).toDouble();

  // Check for deflation.

  while (true) {
    if (N0.value < I0) return;
    if (N0.value == I0) {
      Z[4 * N0.value - 3] = Z[4 * N0.value + PP.value - 3] + SIGMA.value;
      N0.value--;
      continue;
    }
    NN = 4 * N0.value + PP.value;
    if (N0.value != (I0 + 1)) {
      // Check whether E(N0-1) is negligible, 1 eigenvalue.
      if (Z[NN - 5] <= TOL2 * (SIGMA.value + Z[NN - 3]) ||
          Z[NN - 2 * PP.value - 4] <= TOL2 * Z[NN - 7]) {
        Z[4 * N0.value - 3] = Z[4 * N0.value + PP.value - 3] + SIGMA.value;
        N0.value--;
        continue;
      }

      // Check  whether E(N0-2) is negligible, 2 eigenvalues.
      if (Z[NN - 9] > TOL2 * SIGMA.value &&
          Z[NN - 2 * PP.value - 8] > TOL2 * Z[NN - 11]) break;
    }

    if (Z[NN - 3] > Z[NN - 7]) {
      S = Z[NN - 3];
      Z[NN - 3] = Z[NN - 7];
      Z[NN - 7] = S;
    }
    T = HALF * ((Z[NN - 7] - Z[NN - 3]) + Z[NN - 5]);
    if (Z[NN - 5] > Z[NN - 3] * TOL2 && T != ZERO) {
      S = Z[NN - 3] * (Z[NN - 5] / T);
      if (S <= T) {
        S = Z[NN - 3] * (Z[NN - 5] / (T * (ONE + sqrt(ONE + S / T))));
      } else {
        S = Z[NN - 3] * (Z[NN - 5] / (T + sqrt(T) * sqrt(T + S)));
      }
      T = Z[NN - 7] + (S + Z[NN - 5]);
      Z[NN - 3] *= (Z[NN - 7] / T);
      Z[NN - 7] = T;
    }
    Z[4 * N0.value - 7] = Z[NN - 7] + SIGMA.value;
    Z[4 * N0.value - 3] = Z[NN - 3] + SIGMA.value;
    N0.value -= 2;
  }

  if (PP.value == 2) PP.value = 0;

  // Reverse the qd-array, if warranted.

  if (DMIN.value <= ZERO || N0.value < N0IN) {
    if (CBIAS * Z[4 * I0 + PP.value - 3] < Z[4 * N0.value + PP.value - 3]) {
      IPN4 = 4 * (I0 + N0.value);
      for (J4 = 4 * I0; J4 <= 2 * (I0 + N0.value - 1); J4 += 4) {
        TEMP = Z[J4 - 3];
        Z[J4 - 3] = Z[IPN4 - J4 - 3];
        Z[IPN4 - J4 - 3] = TEMP;
        TEMP = Z[J4 - 2];
        Z[J4 - 2] = Z[IPN4 - J4 - 2];
        Z[IPN4 - J4 - 2] = TEMP;
        TEMP = Z[J4 - 1];
        Z[J4 - 1] = Z[IPN4 - J4 - 5];
        Z[IPN4 - J4 - 5] = TEMP;
        TEMP = Z[J4];
        Z[J4] = Z[IPN4 - J4 - 4];
        Z[IPN4 - J4 - 4] = TEMP;
      }
      if (N0.value - I0 <= 4) {
        Z[4 * N0.value + PP.value - 1] = Z[4 * I0 + PP.value - 1];
        Z[4 * N0.value - PP.value] = Z[4 * I0 - PP.value];
      }
      DMIN2.value = min(DMIN2.value, Z[4 * N0.value + PP.value - 1]);
      Z[4 * N0.value + PP.value - 1] = min(Z[4 * N0.value + PP.value - 1],
          min(Z[4 * I0 + PP.value - 1], Z[4 * I0 + PP.value + 3]));
      Z[4 * N0.value - PP.value] = min(Z[4 * N0.value - PP.value],
          min(Z[4 * I0 - PP.value], Z[4 * I0 - PP.value + 4]));
      QMAX.value = max(
          QMAX.value, max(Z[4 * I0 + PP.value - 3], Z[4 * I0 + PP.value + 1]));
      DMIN.value = -ZERO;
    }
  }

  // Choose a shift.

  dlasq4(I0, N0.value, Z, PP.value, N0IN, DMIN.value, DMIN1.value, DMIN2.value,
      DN.value, DN1.value, DN2.value, TAU, TTYPE, G);

  // Call dqds until DMIN.value > 0.

  var success = false;
  while (true) {
    dlasq5(I0, N0.value, Z, PP.value, TAU.value, SIGMA.value, DMIN, DMIN1,
        DMIN2, DN, DN1, DN2, IEEE, EPS);

    NDIV.value += (N0.value - I0 + 2);
    ITER.value++;

    // Check status.

    if (DMIN.value >= ZERO && DMIN1.value >= ZERO) {
      // Success.

      success = true;
      break;
    } else if (DMIN.value < ZERO &&
        DMIN1.value > ZERO &&
        Z[4 * (N0.value - 1) - PP.value] < TOL * (SIGMA.value + DN1.value) &&
        (DN.value).abs() < TOL * SIGMA.value) {
      // Convergence hidden by negative DN.value.

      Z[4 * (N0.value - 1) - PP.value + 2] = ZERO;
      DMIN.value = ZERO;

      success = true;
      break;
    } else if (DMIN.value < ZERO) {
      // TAU.value too big. Select new TAU.value and try again.

      NFAIL.value++;
      if (TTYPE.value < -22) {
        // Failed twice. Play it safe.

        TAU.value = ZERO;
      } else if (DMIN1.value > ZERO) {
        // Late failure. Gives excellent shift.

        TAU.value = (TAU.value + DMIN.value) * (ONE - TWO * EPS);
        TTYPE.value -= 11;
      } else {
        // Early failure. Divide by 4.

        TAU.value = QURTR * TAU.value;
        TTYPE.value -= 12;
      }
      continue;
    } else if (disnan(DMIN.value)) {
      // NaN.

      if (TAU.value == ZERO) {
        break;
      } else {
        TAU.value = ZERO;
        continue;
      }
    } else {
      // Possible underflow. Play it safe.
    }
    break;
  }

  // Risk of underflow.

  if (!success) {
    dlasq6(I0, N0.value, Z, PP.value, DMIN, DMIN1, DMIN2, DN, DN1, DN2);
    NDIV.value += (N0.value - I0 + 2);
    ITER.value++;
    TAU.value = ZERO;
  }

  if (TAU.value < SIGMA.value) {
    DESIG.value += TAU.value;
    T = SIGMA.value + DESIG.value;
    DESIG.value -= (T - SIGMA.value);
  } else {
    T = SIGMA.value + TAU.value;
    DESIG.value = SIGMA.value - (T - TAU.value) + DESIG.value;
  }
  SIGMA.value = T;
}
