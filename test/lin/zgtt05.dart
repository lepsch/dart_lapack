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

  double CABS1(Complex ZDUM) => ZDUM.real.abs() + ZDUM.imaginary.abs();

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
    final XNORM = max(CABS1(X[IMAX][J]), UNFL);
    var DIFF = ZERO;
    for (var I = 1; I <= N; I++) {
      DIFF = max(DIFF, CABS1(X[I][J] - XACT[I][J]));
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
        AXBI = CABS1(B[1][K]) + CABS1(D[1]) * CABS1(X[1][K]);
      } else {
        AXBI = CABS1(B[1][K]) +
            CABS1(D[1]) * CABS1(X[1][K]) +
            CABS1(DU[1]) * CABS1(X[2][K]);
        for (var I = 2; I <= N - 1; I++) {
          final TMP = CABS1(B[I][K]) +
              CABS1(DL[I - 1]) * CABS1(X[I - 1][K]) +
              CABS1(D[I]) * CABS1(X[I][K]) +
              CABS1(DU[I]) * CABS1(X[I + 1][K]);
          AXBI = min(AXBI, TMP);
        }
        final TMP = CABS1(B[N][K]) +
            CABS1(DL[N - 1]) * CABS1(X[N - 1][K]) +
            CABS1(D[N]) * CABS1(X[N][K]);
        AXBI = min(AXBI, TMP);
      }
    } else {
      if (N == 1) {
        AXBI = CABS1(B[1][K]) + CABS1(D[1]) * CABS1(X[1][K]);
      } else {
        AXBI = CABS1(B[1][K]) +
            CABS1(D[1]) * CABS1(X[1][K]) +
            CABS1(DL[1]) * CABS1(X[2][K]);
        for (var I = 2; I <= N - 1; I++) {
          final TMP = CABS1(B[I][K]) +
              CABS1(DU[I - 1]) * CABS1(X[I - 1][K]) +
              CABS1(D[I]) * CABS1(X[I][K]) +
              CABS1(DL[I]) * CABS1(X[I + 1][K]);
          AXBI = min(AXBI, TMP);
        }
        final TMP = CABS1(B[N][K]) +
            CABS1(DU[N - 1]) * CABS1(X[N - 1][K]) +
            CABS1(D[N]) * CABS1(X[N][K]);
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
