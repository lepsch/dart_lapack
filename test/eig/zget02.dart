import 'dart:math';

import 'package:lapack/src/blas/dzasum.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlange.dart';

void zget02(
  final String TRANS,
  final int M,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> X_,
  final int LDX,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
  final A = A_.having(ld: LDA);
  final X = X_.having(ld: LDX);
  final B = B_.having(ld: LDB);
  final RWORK = RWORK_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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
    ANORM = zlange('1', M, N, A, LDA, RWORK);
  } else {
    ANORM = zlange('I', M, N, A, LDA, RWORK);
  }
  if (ANORM <= ZERO) {
    RESID.value = ONE / EPS;
    return;
  }

  // Compute B - op(A)*X and store in B.

  zgemm(TRANS, 'No transpose', N1, NRHS, N2, -Complex.one, A, LDA, X, LDX,
      Complex.one, B, LDB);

  // Compute the maximum over the number of right hand sides of
  //    norm(B - op(A)*X) / ( norm(op(A)) * norm(X) * EPS ) .

  RESID.value = ZERO;
  for (J = 1; J <= NRHS; J++) {
    // 10
    BNORM = dzasum(N1, B(1, J).asArray(), 1);
    XNORM = dzasum(N2, X(1, J).asArray(), 1);
    if (XNORM <= ZERO) {
      RESID.value = ONE / EPS;
    } else {
      RESID.value = max(RESID.value, ((BNORM / ANORM) / XNORM) / EPS);
    }
  } // 10
}
