import 'dart:math';

import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dget04(
  final int N,
  final int NRHS,
  final Matrix<double> X_,
  final int LDX,
  final Matrix<double> XACT_,
  final int LDXACT,
  final double RCOND,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having(ld: LDX);
  final XACT = XACT_.having(ld: LDXACT);
  const ZERO = 0.0;

  // Quick exit if N = 0 or NRHS = 0.

  if (N <= 0 || NRHS <= 0) {
    RESID.value = ZERO;
    return;
  }

  // Exit with RESID.value = 1/EPS if RCOND is invalid.

  final EPS = dlamch('Epsilon');
  if (RCOND < ZERO) {
    RESID.value = 1.0 / EPS;
    return;
  }

  // Compute the maximum of
  //    norm(X - XACT) / ( norm(XACT) * EPS )
  // over all the vectors X and XACT .

  RESID.value = ZERO;
  for (var J = 1; J <= NRHS; J++) {
    final IX = idamax(N, XACT(1, J).asArray(), 1);
    final XNORM = XACT[IX][J].abs();
    var DIFFNM = ZERO;
    for (var I = 1; I <= N; I++) {
      DIFFNM = max(DIFFNM, (X[I][J] - XACT[I][J]).abs());
    }
    if (XNORM <= ZERO) {
      if (DIFFNM > ZERO) RESID.value = 1.0 / EPS;
    } else {
      RESID.value = max(RESID.value, (DIFFNM / XNORM) * RCOND);
    }
  }
  if (RESID.value * EPS < 1.0) RESID.value /= EPS;
}
