import 'dart:math';

import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';

void ztrt05(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
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
  final A = A_.having(ld: LDA);
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
  final UPPER = lsame(UPLO, 'U');
  final NOTRAN = lsame(TRANS, 'N');
  final UNIT = lsame(DIAG, 'U');

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

  // Test 2:  Compute the maximum of BERR / ( (n+1)*EPS + (*) ), where
  // (*) = (n+1)*UNFL / (min_i (abs(A)*abs(X) +abs(b))_i )

  final IFU = UNIT ? 1 : 0;
  for (var K = 1; K <= NRHS; K++) {
    var AXBI = ZERO;
    for (var I = 1; I <= N; I++) {
      var TMP = CABS1(B[I][K]);
      if (UPPER) {
        if (!NOTRAN) {
          for (var J = 1; J <= I - IFU; J++) {
            TMP += CABS1(A[J][I]) * CABS1(X[J][K]);
          }
          if (UNIT) TMP += CABS1(X[I][K]);
        } else {
          if (UNIT) TMP += CABS1(X[I][K]);
          for (var J = I + IFU; J <= N; J++) {
            TMP += CABS1(A[I][J]) * CABS1(X[J][K]);
          }
        }
      } else {
        if (NOTRAN) {
          for (var J = 1; J <= I - IFU; J++) {
            TMP += CABS1(A[I][J]) * CABS1(X[J][K]);
          }
          if (UNIT) TMP += CABS1(X[I][K]);
        } else {
          if (UNIT) TMP += CABS1(X[I][K]);
          for (var J = I + IFU; J <= N; J++) {
            TMP += CABS1(A[J][I]) * CABS1(X[J][K]);
          }
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
