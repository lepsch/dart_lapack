import 'package:lapack/src/blas/zdotu.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlansp.dart';

void zspt03(
  final String UPLO,
  final int N,
  final Array<Complex> A_,
  final Array<Complex> AINV_,
  final Matrix<Complex> WORK_,
  final int LDW,
  final Array<double> RWORK_,
  final Box<double> RCOND,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having();
  final AINV = AINV_.having();
  final WORK = WORK_.having(ld: LDW);
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if N = 0.

  if (N <= 0) {
    RCOND.value = ONE;
    RESID.value = ZERO;
    return;
  }

  // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

  final EPS = dlamch('Epsilon');
  final ANORM = zlansp('1', UPLO, N, A, RWORK);
  final AINVNM = zlansp('1', UPLO, N, AINV, RWORK);
  if (ANORM <= ZERO || AINVNM <= ZERO) {
    RCOND.value = ZERO;
    RESID.value = ONE / EPS;
    return;
  }
  RCOND.value = (ONE / ANORM) / AINVNM;

  // Case where both A and AINV are upper triangular:
  // Each element of - A * AINV is computed by taking the dot product
  // of a row of A with a column of AINV.

  if (lsame(UPLO, 'U')) {
    for (var I = 1; I <= N; I++) {
      final ICOL = ((I - 1) * I) ~/ 2 + 1;

      // Code when J <= I

      for (var J = 1; J <= I; J++) {
        var JCOL = ((J - 1) * J) ~/ 2 + 1;
        var T = zdotu(J, A(ICOL), 1, AINV(JCOL), 1);
        JCOL += 2 * J - 1;
        var KCOL = ICOL - 1;
        for (var K = J + 1; K <= I; K++) {
          T += A[KCOL + K] * AINV[JCOL];
          JCOL += K;
        }
        KCOL += 2 * I;
        for (var K = I + 1; K <= N; K++) {
          T += A[KCOL] * AINV[JCOL];
          KCOL += K;
          JCOL += K;
        }
        WORK[I][J] = -T;
      }

      // Code when J > I

      for (var J = I + 1; J <= N; J++) {
        var JCOL = ((J - 1) * J) ~/ 2 + 1;
        var T = zdotu(I, A(ICOL), 1, AINV(JCOL), 1);
        JCOL--;
        var KCOL = ICOL + 2 * I - 1;
        for (var K = I + 1; K <= J; K++) {
          T += A[KCOL] * AINV[JCOL + K];
          KCOL += K;
        }
        JCOL += 2 * J;
        for (var K = J + 1; K <= N; K++) {
          T += A[KCOL] * AINV[JCOL];
          KCOL += K;
          JCOL += K;
        }
        WORK[I][J] = -T;
      }
    }
  } else {
    // Case where both A and AINV are lower triangular

    final NALL = (N * (N + 1)) ~/ 2;
    for (var I = 1; I <= N; I++) {
      // Code when J <= I

      var ICOL = NALL - ((N - I + 1) * (N - I + 2)) ~/ 2 + 1;
      for (var J = 1; J <= I; J++) {
        var JCOL = NALL - ((N - J) * (N - J + 1)) ~/ 2 - (N - I);
        var T = zdotu(N - I + 1, A(ICOL), 1, AINV(JCOL), 1);
        var KCOL = I;
        JCOL = J;
        for (var K = 1; K <= J - 1; K++) {
          T += A[KCOL] * AINV[JCOL];
          JCOL += N - K;
          KCOL += N - K;
        }
        JCOL -= J;
        for (var K = J; K <= I - 1; K++) {
          T += A[KCOL] * AINV[JCOL + K];
          KCOL += N - K;
        }
        WORK[I][J] = -T;
      }

      // Code when J > I

      ICOL = NALL - ((N - I) * (N - I + 1)) ~/ 2;
      for (var J = I + 1; J <= N; J++) {
        var JCOL = NALL - ((N - J + 1) * (N - J + 2)) ~/ 2 + 1;
        var T = zdotu(N - J + 1, A(ICOL - N + J), 1, AINV(JCOL), 1);
        var KCOL = I;
        JCOL = J;
        for (var K = 1; K <= I - 1; K++) {
          T += A[KCOL] * AINV[JCOL];
          JCOL += N - K;
          KCOL += N - K;
        }
        KCOL -= I;
        for (var K = I; K <= J - 1; K++) {
          T += A[KCOL + K] * AINV[JCOL];
          JCOL += N - K;
        }
        WORK[I][J] = -T;
      }
    }
  }

  // Add the identity matrix to WORK .

  for (var I = 1; I <= N; I++) {
    WORK[I][I] += Complex.one;
  }

  // Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)

  RESID.value = zlange('1', N, N, WORK, LDW, RWORK);

  RESID.value = ((RESID.value * RCOND.value) / EPS) / N;
}
