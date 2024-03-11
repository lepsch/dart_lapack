import 'dart:math';

import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';

void dppt05(
  final String UPLO,
  final int N,
  final int NRHS,
  final Array<double> AP_,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> X_,
  final int LDX,
  final Matrix<double> XACT_,
  final int LDXACT,
  final Array<double> FERR_,
  final Array<double> BERR_,
  final Array<double> RESLTS_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final XACT = XACT_.having(ld: LDXACT);
  final AP = AP_.having();
  final FERR = FERR_.having();
  final BERR = BERR_.having();
  final RESLTS = RESLTS_.having();

  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if N = 0 or NRHS = 0.

  if (N <= 0 || NRHS <= 0) {
    RESLTS[1] = ZERO;
    RESLTS[2] = ZERO;
    return;
  }

  final EPS = dlamch('Epsilon');
  final UNFL = dlamch('Safe minimum');
  final OVFL = ONE / UNFL;
  final UPPER = lsame(UPLO, 'U');

  // Test 1:  Compute the maximum of
  //    norm(X - XACT) / ( norm(X) * FERR )
  // over all the vectors X and XACT using the infinity-norm.

  var ERRBND = ZERO;
  for (var J = 1; J <= NRHS; J++) {
    final IMAX = idamax(N, X(1, J).asArray(), 1);
    final XNORM = max(X[IMAX][J].abs(), UNFL);
    var DIFF = ZERO;
    for (var I = 1; I <= N; I++) {
      DIFF = max(DIFF, (X[I][J] - XACT[I][J]).abs());
    }

    if (XNORM > ONE) {
      //
    } else if (DIFF <= OVFL * XNORM) {
      //
    } else {
      ERRBND = ONE / EPS;
      continue;
    }

    if (DIFF / XNORM <= FERR[J]) {
      ERRBND = max(ERRBND, (DIFF / XNORM) / FERR[J]);
    } else {
      ERRBND = ONE / EPS;
    }
  }
  RESLTS[1] = ERRBND;

  // Test 2:  Compute the maximum of BERR / ( (n+1)*EPS + (*) ), where
  // (*) = (n+1)*UNFL / (min_i (abs(A)*abs(X) +abs(b))_i )

  for (var K = 1; K <= NRHS; K++) {
    double AXBI = 0;
    for (var I = 1; I <= N; I++) {
      var TMP = B[I][K].abs();
      if (UPPER) {
        var JC = ((I - 1) * I) ~/ 2;
        for (var J = 1; J <= I; J++) {
          TMP += (AP[JC + J]).abs() * (X[J][K]).abs();
        }
        JC = JC + I;
        for (var J = I + 1; J <= N; J++) {
          TMP += (AP[JC]).abs() * (X[J][K]).abs();
          JC = JC + J;
        }
      } else {
        var JC = I;
        for (var J = 1; J <= I - 1; J++) {
          TMP += (AP[JC]).abs() * (X[J][K]).abs();
          JC = JC + N - J;
        }
        for (var J = I; J <= N; J++) {
          TMP += (AP[JC + J - I]).abs() * (X[J][K]).abs();
        }
      }
      if (I == 1) {
        AXBI = TMP;
      } else {
        AXBI = min(AXBI, TMP);
      }
    }
    final TMP =
        BERR[K] / ((N + 1) * EPS + (N + 1) * UNFL / max(AXBI, (N + 1) * UNFL));
    if (K == 1) {
      RESLTS[2] = TMP;
    } else {
      RESLTS[2] = max(RESLTS[2], TMP);
    }
  }
}
