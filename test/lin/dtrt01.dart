import 'package:lapack/src/blas/dtrmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlantr.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';

void dtrt01(
  final String UPLO,
  final String DIAG,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> AINV_,
  final int LDAINV,
  final Box<double> RCOND,
  final Array<double> WORK_,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final AINV = AINV_.having(ld: LDAINV);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if N = 0

  if (N <= 0) {
    RCOND.value = ONE;
    RESID.value = ZERO;
    return;
  }

  // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

  final EPS = dlamch('Epsilon');
  final ANORM = dlantr('1', UPLO, DIAG, N, N, A, LDA, WORK);
  final AINVNM = dlantr('1', UPLO, DIAG, N, N, AINV, LDAINV, WORK);
  if (ANORM <= ZERO || AINVNM <= ZERO) {
    RCOND.value = ZERO;
    RESID.value = ONE / EPS;
    return;
  }
  RCOND.value = (ONE / ANORM) / AINVNM;

  // Set the diagonal of AINV to 1 if AINV has unit diagonal.

  if (lsame(DIAG, 'U')) {
    for (var J = 1; J <= N; J++) {
      AINV[J][J] = ONE;
    }
  }

  // Compute A * AINV, overwriting AINV.

  if (lsame(UPLO, 'U')) {
    for (var J = 1; J <= N; J++) {
      dtrmv('Upper', 'No transpose', DIAG, J, A, LDA, AINV(1, J).asArray(), 1);
    }
  } else {
    for (var J = 1; J <= N; J++) {
      dtrmv('Lower', 'No transpose', DIAG, N - J + 1, A(J, J), LDA,
          AINV(J, J).asArray(), 1);
    }
  }

  // Subtract 1 from each diagonal element to form A*AINV - I.

  for (var J = 1; J <= N; J++) {
    AINV[J][J] -= ONE;
  }

  // Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)

  RESID.value = dlantr('1', UPLO, 'Non-unit', N, N, AINV, LDAINV, WORK);

  RESID.value = ((RESID.value * RCOND.value) / N) / EPS;
}
