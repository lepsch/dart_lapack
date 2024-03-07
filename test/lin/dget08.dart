import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dget08(
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
  final A = A_.having(ld: LDA);
  final X = X_.having(ld: LDX);
  final B = B_.having(ld: LDB);
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if M = 0 or N = 0 or NRHS = 0

  if (M <= 0 || N <= 0 || NRHS == 0) {
    RESID.value = ZERO;
    return;
  }

  final (N1, N2) = lsame(TRANS, 'T') || lsame(TRANS, 'C') ? (N, M) : (M, N);

  // Exit with RESID.value = 1/EPS if ANORM = 0.

  final EPS = dlamch('Epsilon');
  final ANORM = dlange('I', N1, N2, A, LDA, RWORK);
  if (ANORM <= ZERO) {
    RESID.value = ONE / EPS;
    return;
  }

  // Compute  B - A*X  (or  B - A'*X ) and store in B.

  dgemm(TRANS, 'No transpose', N1, NRHS, N2, -ONE, A, LDA, X, LDX, ONE, B, LDB);

  // Compute the maximum over the number of right hand sides of
  //    norm(B - A*X) / ( norm(A) * norm(X) * EPS ) .

  RESID.value = ZERO;
  for (var J = 1; J <= NRHS; J++) {
    final BNORM = B[idamax(N1, B(1, J).asArray(), 1)][J].abs();
    final XNORM = X[idamax(N2, X(1, J).asArray(), 1)][J].abs();
    if (XNORM <= ZERO) {
      RESID.value = ONE / EPS;
    } else {
      RESID.value = max(RESID.value, ((BNORM / ANORM) / XNORM) / EPS);
    }
  }
}
