import 'dart:math';

import 'package:lapack/src/blas/dasum.dart';
import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dtrmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlantr.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';

void dtrt02(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> X_,
  final int LDX,
  final Matrix<double> B_,
  final int LDB,
  final Array<double> WORK_,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final X = X_.having(ld: LDX);
  final B = B_.having(ld: LDB);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if N = 0 or NRHS = 0

  if (N <= 0 || NRHS <= 0) {
    RESID.value = ZERO;
    return;
  }

  // Compute the 1-norm of op(A).

  final double ANORM;
  if (lsame(TRANS, 'N')) {
    ANORM = dlantr('1', UPLO, DIAG, N, N, A, LDA, WORK);
  } else {
    ANORM = dlantr('I', UPLO, DIAG, N, N, A, LDA, WORK);
  }

  // Exit with RESID.value = 1/EPS if ANORM = 0.

  final EPS = dlamch('Epsilon');
  if (ANORM <= ZERO) {
    RESID.value = ONE / EPS;
    return;
  }

  // Compute the maximum over the number of right hand sides of
  //    norm(op(A)*X - B) / ( norm(op(A)) * norm(X) * EPS )

  RESID.value = ZERO;
  for (var J = 1; J <= NRHS; J++) {
    dcopy(N, X(1, J).asArray(), 1, WORK, 1);
    dtrmv(UPLO, TRANS, DIAG, N, A, LDA, WORK, 1);
    daxpy(N, -ONE, B(1, J).asArray(), 1, WORK, 1);
    final BNORM = dasum(N, WORK, 1);
    final XNORM = dasum(N, X(1, J).asArray(), 1);
    if (XNORM <= ZERO) {
      RESID.value = ONE / EPS;
    } else {
      RESID.value = max(RESID.value, ((BNORM / ANORM) / XNORM) / EPS);
    }
  }
}
