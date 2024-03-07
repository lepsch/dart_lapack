import 'dart:math';

import 'package:lapack/src/blas/dasum.dart';
import 'package:lapack/src/blas/dgbmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';

void dgbt02(
  final String TRANS,
  final int M,
  final int N,
  final int KL,
  final int KU,
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

  // Quick return if N = 0 pr NRHS = 0

  if (M <= 0 || N <= 0 || NRHS <= 0) {
    RESID.value = ZERO;
    return;
  }

  // Exit with RESID = 1/EPS if ANORM = 0.

  final EPS = dlamch('Epsilon');
  var ANORM = ZERO;
  if (lsame(TRANS, 'N')) {
    // Find norm1(A).

    final KD = KU + 1;
    for (var J = 1; J <= N; J++) {
      final I1 = max(KD + 1 - J, 1);
      final I2 = min(KD + M - J, KL + KD);
      if (I2 >= I1) {
        final TEMP = dasum(I2 - I1 + 1, A(I1, J).asArray(), 1);
        if (ANORM < TEMP || disnan(TEMP)) ANORM = TEMP;
      }
    }
  } else {
    // Find normI(A).

    for (var I1 = 1; I1 <= M; I1++) {
      RWORK[I1] = ZERO;
    }
    for (var J = 1; J <= N; J++) {
      final KD = KU + 1 - J;
      for (var I1 = max(1, J - KU); I1 <= min(M, J + KL); I1++) {
        RWORK[I1] = RWORK[I1] + (A[KD + I1][J]).abs();
      }
    }
    for (var I1 = 1; I1 <= M; I1++) {
      final TEMP = RWORK[I1];
      if (ANORM < TEMP || disnan(TEMP)) ANORM = TEMP;
    }
  }
  if (ANORM <= ZERO) {
    RESID.value = ONE / EPS;
    return;
  }

  final N1 = lsame(TRANS, 'T') || lsame(TRANS, 'C') ? N : M;

  // Compute B - op(A)*X

  for (var J = 1; J <= NRHS; J++) {
    dgbmv(TRANS, M, N, KL, KU, -ONE, A, LDA, X(1, J).asArray(), 1, ONE,
        B(1, J).asArray(), 1);
  }

  // Compute the maximum over the number of right hand sides of
  //    norm(B - op(A)*X) / ( norm(op(A)) * norm(X) * EPS ).

  RESID.value = ZERO;
  for (var J = 1; J <= NRHS; J++) {
    final BNORM = dasum(N1, B(1, J).asArray(), 1);
    final XNORM = dasum(N1, X(1, J).asArray(), 1);
    if (XNORM <= ZERO) {
      RESID.value = ONE / EPS;
    } else {
      RESID.value = max(RESID.value, ((BNORM / ANORM) / XNORM) / EPS);
    }
  }
}
