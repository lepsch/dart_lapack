import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zhpmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlanhp.dart';

void zppt03(
  final String UPLO,
  final int N,
  final Array<Complex> A_,
  final Array<Complex> AINV_,
  final Matrix<Complex> WORK_,
  final int LDWORK,
  final Array<double> RWORK_,
  final Box<double> RCOND,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having();
  final AINV = AINV_.having();
  final WORK = WORK_.having(ld: LDWORK);
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
  final ANORM = zlanhp('1', UPLO, N, A, RWORK);
  final AINVNM = zlanhp('1', UPLO, N, AINV, RWORK);
  if (ANORM <= ZERO || AINVNM <= ZERO) {
    RCOND.value = ZERO;
    RESID.value = ONE / EPS;
    return;
  }
  RCOND.value = (ONE / ANORM) / AINVNM;

  // UPLO = 'U':
  // Copy the leading N-1 x N-1 submatrix of AINV to WORK(1:N,2:N) and
  // expand it to a full matrix, then multiply by A one column at a
  // time, moving the result one column to the left.

  if (lsame(UPLO, 'U')) {
    // Copy AINV

    var JJ = 1;
    for (var J = 1; J <= N - 1; J++) {
      zcopy(J, AINV(JJ), 1, WORK(1, J + 1).asArray(), 1);
      for (var I = 1; I <= J - 1; I++) {
        WORK[J][I + 1] = AINV[JJ + I - 1].conjugate();
      }
      JJ += J;
    }
    JJ = ((N - 1) * N) ~/ 2 + 1;
    for (var I = 1; I <= N - 1; I++) {
      WORK[N][I + 1] = AINV[JJ + I - 1].conjugate();
    }

    // Multiply by A

    for (var J = 1; J <= N - 1; J++) {
      zhpmv('Upper', N, -Complex.one, A, WORK(1, J + 1).asArray(), 1,
          Complex.zero, WORK(1, J).asArray(), 1);
    }
    zhpmv('Upper', N, -Complex.one, A, AINV(JJ), 1, Complex.zero,
        WORK(1, N).asArray(), 1);

    // UPLO = 'L':
    // Copy the trailing N-1 x N-1 submatrix of AINV to WORK(1:N,1:N-1)
    // and multiply by A, moving each column to the right.
  } else {
    // Copy AINV

    for (var I = 1; I <= N - 1; I++) {
      WORK[1][I] = AINV[I + 1].conjugate();
    }
    var JJ = N + 1;
    for (var J = 2; J <= N; J++) {
      zcopy(N - J + 1, AINV(JJ), 1, WORK(J, J - 1).asArray(), 1);
      for (var I = 1; I <= N - J; I++) {
        WORK[J][J + I - 1] = AINV[JJ + I].conjugate();
      }
      JJ += N - J + 1;
    }

    // Multiply by A

    for (var J = N; J >= 2; J--) {
      zhpmv('Lower', N, -Complex.one, A, WORK(1, J - 1).asArray(), 1,
          Complex.zero, WORK(1, J).asArray(), 1);
    }
    zhpmv('Lower', N, -Complex.one, A, AINV(1), 1, Complex.zero,
        WORK(1, 1).asArray(), 1);
  }

  // Add the identity matrix to WORK .

  for (var I = 1; I <= N; I++) {
    WORK[I][I] += Complex.one;
  }

  // Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)

  RESID.value = zlange('1', N, N, WORK, LDWORK, RWORK);

  RESID.value = ((RESID.value * RCOND.value) / EPS) / N;
}
