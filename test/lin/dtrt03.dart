import 'dart:math';

import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dtrmv.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';

void dtrt03(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final double SCALE,
  final Array<double> CNORM_,
  final double TSCAL,
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
  final CNORM = CNORM_.having();
  final WORK = WORK_.having();
  const ONE = 1.0, ZERO = 0.0;

  // Quick exit if N = 0

  if (N <= 0 || NRHS <= 0) {
    RESID.value = ZERO;
    return;
  }
  final EPS = dlamch('Epsilon');
  final SMLNUM = dlamch('Safe minimum');

  // Compute the norm of the triangular matrix A using the column
  // norms already computed by DLATRS.

  var TNORM = ZERO;
  if (lsame(DIAG, 'N')) {
    for (var J = 1; J <= N; J++) {
      TNORM = max(TNORM, TSCAL * (A[J][J]).abs() + CNORM[J]);
    }
  } else {
    for (var J = 1; J <= N; J++) {
      TNORM = max(TNORM, TSCAL + CNORM[J]);
    }
  }

  // Compute the maximum over the number of right hand sides of
  //    norm(op(A)*x - s*b) / ( norm(op(A)) * norm(x) * EPS ).

  RESID.value = ZERO;
  for (var J = 1; J <= NRHS; J++) {
    dcopy(N, X(1, J).asArray(), 1, WORK, 1);
    var IX = idamax(N, WORK, 1);
    var XNORM = max(ONE, (X[IX][J]).abs());
    final XSCAL = (ONE / XNORM) / N;
    dscal(N, XSCAL, WORK, 1);
    dtrmv(UPLO, TRANS, DIAG, N, A, LDA, WORK, 1);
    daxpy(N, -SCALE * XSCAL, B(1, J).asArray(), 1, WORK, 1);
    IX = idamax(N, WORK, 1);
    var ERR = TSCAL * (WORK[IX]).abs();
    IX = idamax(N, X(1, J).asArray(), 1);
    XNORM = (X[IX][J]).abs();
    if (ERR * SMLNUM <= XNORM) {
      if (XNORM > ZERO) ERR = ERR / XNORM;
    } else {
      if (ERR > ZERO) ERR = ONE / EPS;
    }
    if (ERR * SMLNUM <= TNORM) {
      if (TNORM > ZERO) ERR = ERR / TNORM;
    } else {
      if (ERR > ZERO) ERR = ONE / EPS;
    }
    RESID.value = max(RESID.value, ERR);
  }
}
