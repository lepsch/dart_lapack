import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dsycon.dart';
import 'package:lapack/src/dsyrfs.dart';
import 'package:lapack/src/dsytrf.dart';
import 'package:lapack/src/dsytrs.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dsysvx(
  final String FACT,
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> AF_,
  final int LDAF,
  final Array<int> IPIV_,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> X_,
  final int LDX,
  final Box<double> RCOND,
  final Array<double> FERR_,
  final Array<double> BERR_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDAF);
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final FERR = FERR_.having();
  final BERR = BERR_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0;
  bool LQUERY, NOFACT;
  int LWKMIN, LWKOPT = 0, NB;
  double ANORM;

  // Test the input parameters.

  INFO.value = 0;
  NOFACT = lsame(FACT, 'N');
  LQUERY = (LWORK == -1);
  LWKMIN = max(1, 3 * N);
  if (!NOFACT && !lsame(FACT, 'F')) {
    INFO.value = -1;
  } else if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (NRHS < 0) {
    INFO.value = -4;
  } else if (LDA < max(1, N)) {
    INFO.value = -6;
  } else if (LDAF < max(1, N)) {
    INFO.value = -8;
  } else if (LDB < max(1, N)) {
    INFO.value = -11;
  } else if (LDX < max(1, N)) {
    INFO.value = -13;
  } else if (LWORK < LWKMIN && !LQUERY) {
    INFO.value = -18;
  }

  if (INFO.value == 0) {
    LWKOPT = LWKMIN;
    if (NOFACT) {
      NB = ilaenv(1, 'DSYTRF', UPLO, N, -1, -1, -1);
      LWKOPT = max(LWKOPT, N * NB);
    }
    WORK[1] = LWKOPT.toDouble();
  }

  if (INFO.value != 0) {
    xerbla('DSYSVX', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  if (NOFACT) {
    // Compute the factorization A = U*D*U**T or A = L*D*L**T.

    dlacpy(UPLO, N, N, A, LDA, AF, LDAF);
    dsytrf(UPLO, N, AF, LDAF, IPIV, WORK, LWORK, INFO);

    // Return if INFO.value is non-zero.

    if (INFO.value > 0) {
      RCOND.value = ZERO;
      return;
    }
  }

  // Compute the norm of the matrix A.

  ANORM = dlansy('I', UPLO, N, A, LDA, WORK);

  // Compute the reciprocal of the condition number of A.

  dsycon(UPLO, N, AF, LDAF, IPIV, ANORM, RCOND, WORK, IWORK, INFO);

  // Compute the solution vectors X.

  dlacpy('Full', N, NRHS, B, LDB, X, LDX);
  dsytrs(UPLO, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO);

  // Use iterative refinement to improve the computed solutions and
  // compute error bounds and backward error estimates for them.

  dsyrfs(UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR,
      WORK, IWORK, INFO);

  // Set INFO.value = N+1 if the matrix is singular to working precision.

  if (RCOND.value < dlamch('Epsilon')) INFO.value = N + 1;

  WORK[1] = LWKOPT.toDouble();
}
