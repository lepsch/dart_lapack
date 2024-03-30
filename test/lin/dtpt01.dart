import 'package:lapack/src/blas/dtpmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlantp.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';

void dtpt01(
  final String UPLO,
  final String DIAG,
  final int N,
  final Array<double> AP_,
  final Array<double> AINVP_,
  final Box<double> RCOND,
  final Array<double> WORK_,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.having();
  final AINVP = AINVP_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if N = 0.

  if (N <= 0) {
    RCOND.value = ONE;
    RESID.value = ZERO;
    return;
  }

  // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

  final EPS = dlamch('Epsilon');
  final ANORM = dlantp('1', UPLO, DIAG, N, AP, WORK);
  final AINVNM = dlantp('1', UPLO, DIAG, N, AINVP, WORK);
  if (ANORM <= ZERO || AINVNM <= ZERO) {
    RCOND.value = ZERO;
    RESID.value = ONE / EPS;
    return;
  }
  RCOND.value = (ONE / ANORM) / AINVNM;

  // Compute A * AINV, overwriting AINV.

  final UNITD = lsame(DIAG, 'U');
  if (lsame(UPLO, 'U')) {
    var JC = 1;
    for (var J = 1; J <= N; J++) {
      if (UNITD) AINVP[JC + J - 1] = ONE;

      // Form the j-th column of A*AINV

      dtpmv('Upper', 'No transpose', DIAG, J, AP, AINVP(JC), 1);

      // Subtract 1 from the diagonal

      AINVP[JC + J - 1] -= ONE;
      JC += J;
    }
  } else {
    var JC = 1;
    for (var J = 1; J <= N; J++) {
      if (UNITD) AINVP[JC] = ONE;

      // Form the j-th column of A*AINV

      dtpmv('Lower', 'No transpose', DIAG, N - J + 1, AP(JC), AINVP(JC), 1);

      // Subtract 1 from the diagonal

      AINVP[JC] -= ONE;
      JC += N - J + 1;
    }
  }

  // Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)

  RESID.value = dlantp('1', UPLO, 'Non-unit', N, AINVP, WORK);

  RESID.value = ((RESID.value * RCOND.value) / N) / EPS;
}
