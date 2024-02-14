import 'dart:math';

import 'package:lapack/src/blas/dasum.dart';
import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dget02(
  final String TRANS,
  final int M,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> X_,
  final int LDX,
  final Matrix<double> B_,
  final int LDB,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final X = X_.dim(LDX);
  final B = B_.dim(LDB);
  final RWORK = RWORK_.dim();
  const ZERO = 0.0, ONE = 1.0;
  int J, N1, N2;
  double ANORM, BNORM, EPS, XNORM;

  // Quick exit if M = 0 or N = 0 or NRHS = 0

  if (M <= 0 || N <= 0 || NRHS == 0) {
    RESID.value = ZERO;
    return;
  }

  if (lsame(TRANS, 'T') || lsame(TRANS, 'C')) {
    N1 = N;
    N2 = M;
  } else {
    N1 = M;
    N2 = N;
  }

  // Exit with RESID.value = 1/EPS if ANORM = 0.

  EPS = dlamch('Epsilon');
  if (lsame(TRANS, 'N')) {
    ANORM = dlange('1', M, N, A, LDA, RWORK);
  } else {
    ANORM = dlange('I', M, N, A, LDA, RWORK);
  }
  if (ANORM <= ZERO) {
    RESID.value = ONE / EPS;
    return;
  }

  // Compute B - op(A)*X and store in B.

  dgemm(TRANS, 'No transpose', N1, NRHS, N2, -ONE, A, LDA, X, LDX, ONE, B, LDB);

  // Compute the maximum over the number of right hand sides of
  //    norm(B - op(A)*X) / ( norm(op(A)) * norm(X) * EPS ) .

  RESID.value = ZERO;
  for (J = 1; J <= NRHS; J++) {
    BNORM = dasum(N1, B(1, J).asArray(), 1);
    XNORM = dasum(N2, X(1, J).asArray(), 1);
    if (XNORM <= ZERO) {
      RESID.value = ONE / EPS;
    } else {
      RESID.value = max(RESID.value, ((BNORM / ANORM) / XNORM) / EPS);
    }
  }
}
