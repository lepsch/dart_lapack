import 'dart:math';

import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';

void zgtt05(
  final String TRANS,
  final int N,
  final int NRHS,
  final Array<Complex> DL_,
  final Array<Complex> D_,
  final Array<Complex> DU_,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> X_,
  final int LDX,
  final Matrix<Complex> XACT_,
  final int LDXACT,
  final Array<double> FERR_,
  final Array<double> BERR_,
  final Array<double> RESLTS_,
) {
  // -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DL = DL_.having();
  final D = D_.having();
  final DU = DU_.having();
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final XACT = XACT_.having(ld: LDXACT);
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
  final NOTRAN = lsame(TRANS, 'N');
  final NZ = 4;

  // Test 1:  Compute the maximum of
  //    norm(X - XACT) / ( norm(X) * FERR )
  // over all the vectors X and XACT using the infinity-norm.

  var ERRBND = ZERO;
  for (var J = 1; J <= NRHS; J++) {
    final IMAX = izamax(N, X(1, J).asArray(), 1);
    final XNORM = max(X[IMAX][J].cabs1(), UNFL);
    var DIFF = ZERO;
    for (var I = 1; I <= N; I++) {
      DIFF = max(DIFF, (X[I][J] - XACT[I][J]).cabs1());
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

  // Test 2:  Compute the maximum of BERR / ( NZ*EPS + (*) ), where
  // (*) = NZ*UNFL / (min_i (abs(op(A))*abs(X) +abs(b))_i )

  for (var K = 1; K <= NRHS; K++) {
    double AXBI;
    if (NOTRAN) {
      if (N == 1) {
        AXBI = B[1][K].cabs1() + D[1].cabs1() * X[1][K].cabs1();
      } else {
        AXBI = B[1][K].cabs1() +
            D[1].cabs1() * X[1][K].cabs1() +
            DU[1].cabs1() * X[2][K].cabs1();
        for (var I = 2; I <= N - 1; I++) {
          final TMP = B[I][K].cabs1() +
              DL[I - 1].cabs1() * X[I - 1][K].cabs1() +
              D[I].cabs1() * X[I][K].cabs1() +
              DU[I].cabs1() * X[I + 1][K].cabs1();
          AXBI = min(AXBI, TMP);
        }
        final TMP = B[N][K].cabs1() +
            DL[N - 1].cabs1() * X[N - 1][K].cabs1() +
            D[N].cabs1() * X[N][K].cabs1();
        AXBI = min(AXBI, TMP);
      }
    } else {
      if (N == 1) {
        AXBI = B[1][K].cabs1() + D[1].cabs1() * X[1][K].cabs1();
      } else {
        AXBI = B[1][K].cabs1() +
            D[1].cabs1() * X[1][K].cabs1() +
            DL[1].cabs1() * X[2][K].cabs1();
        for (var I = 2; I <= N - 1; I++) {
          final TMP = B[I][K].cabs1() +
              DU[I - 1].cabs1() * X[I - 1][K].cabs1() +
              D[I].cabs1() * X[I][K].cabs1() +
              DL[I].cabs1() * X[I + 1][K].cabs1();
          AXBI = min(AXBI, TMP);
        }
        final TMP = B[N][K].cabs1() +
            DU[N - 1].cabs1() * X[N - 1][K].cabs1() +
            D[N].cabs1() * X[N][K].cabs1();
        AXBI = min(AXBI, TMP);
      }
    }
    final TMP = BERR[K] / (NZ * EPS + NZ * UNFL / max(AXBI, NZ * UNFL));
    if (K == 1) {
      RESLTS[2] = TMP;
    } else {
      RESLTS[2] = max(RESLTS[2], TMP);
    }
  }
}
