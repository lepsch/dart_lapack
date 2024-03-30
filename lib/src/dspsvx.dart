import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlansp.dart';
import 'package:lapack/src/dspcon.dart';
import 'package:lapack/src/dsprfs.dart';
import 'package:lapack/src/dsptrf.dart';
import 'package:lapack/src/dsptrs.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dspsvx(
  final String FACT,
  final String UPLO,
  final int N,
  final int NRHS,
  final Array<double> AP_,
  final Array<double> AFP_,
  final Array<int> IPIV_,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> X_,
  final int LDX,
  final Box<double> RCOND,
  final Array<double> FERR_,
  final Array<double> BERR_,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final AP = AP_.having();
  final AFP = AFP_.having();
  final FERR = FERR_.having();
  final BERR = BERR_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0;
  bool NOFACT;
  double ANORM;

  // Test the input parameters.

  INFO.value = 0;
  NOFACT = lsame(FACT, 'N');
  if (!NOFACT && !lsame(FACT, 'F')) {
    INFO.value = -1;
  } else if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (NRHS < 0) {
    INFO.value = -4;
  } else if (LDB < max(1, N)) {
    INFO.value = -9;
  } else if (LDX < max(1, N)) {
    INFO.value = -11;
  }
  if (INFO.value != 0) {
    xerbla('DSPSVX', -INFO.value);
    return;
  }

  if (NOFACT) {
    // Compute the factorization A = U*D*U**T or A = L*D*L**T.

    dcopy(N * (N + 1) ~/ 2, AP, 1, AFP, 1);
    dsptrf(UPLO, N, AFP, IPIV, INFO);

    // Return if INFO is non-zero.

    if (INFO.value > 0) {
      RCOND.value = ZERO;
      return;
    }
  }

  // Compute the norm of the matrix A.

  ANORM = dlansp('I', UPLO, N, AP, WORK);

  // Compute the reciprocal of the condition number of A.

  dspcon(UPLO, N, AFP, IPIV, ANORM, RCOND, WORK, IWORK, INFO);

  // Compute the solution vectors X.

  dlacpy('Full', N, NRHS, B, LDB, X, LDX);
  dsptrs(UPLO, N, NRHS, AFP, IPIV, X, LDX, INFO);

  // Use iterative refinement to improve the computed solutions and
  // compute error bounds and backward error estimates for them.

  dsprfs(UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK,
      INFO);

  // Set INFO = N+1 if the matrix is singular to working precision.

  if (RCOND.value < dlamch('Epsilon')) INFO.value = N + 1;
}
