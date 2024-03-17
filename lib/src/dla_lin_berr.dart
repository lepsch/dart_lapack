import 'dart:math';

import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dla_lin_berr(
  final int N,
  final int NZ,
  final int NRHS,
  final Matrix<double> RES_,
  final Matrix<double> AYB_,
  final Array<double> BERR_,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final RES = RES_.having(ld: N);
  final AYB = AYB_.having(ld: N);
  final BERR = BERR_.having(length: NRHS);
  double TMP;
  int I, J;

  // Adding SAFE1 to the numerator guards against spuriously zero
  // residuals.  A similar safeguard is in the SLA_yyAMV routine used
  // to compute AYB.

  var SAFE1 = dlamch('Safe minimum');
  SAFE1 = (NZ + 1) * SAFE1;

  for (J = 1; J <= NRHS; J++) {
    BERR[J] = 0.0;
    for (I = 1; I <= N; I++) {
      if (AYB[I][J] != 0.0) {
        TMP = (SAFE1 + RES[I][J].abs()) / AYB[I][J];
        BERR[J] = max(BERR[J], TMP);
      }

      // If AYB is exactly 0.0 (and if computed by SLA_yyAMV), then we know
      // the true residual also must be exactly 0.0.
    }
  }
}
