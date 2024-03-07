import 'dart:math';

import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dtrmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlaswp.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dget01(
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> AFAC_,
  final int LDAFAC,
  final Array<int> IPIV_,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final AFAC = AFAC_.having(ld: LDAFAC);
  final IPIV = IPIV_.having();
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if M = 0 or N = 0.

  if (M <= 0 || N <= 0) {
    RESID.value = ZERO;
    return;
  }

  // Determine EPS and the norm of A.

  final EPS = dlamch('Epsilon');
  final ANORM = dlange('1', M, N, A, LDA, RWORK);

  // Compute the product L*U and overwrite AFAC with the result.
  // A column at a time of the product is obtained, starting with
  // column N.

  for (var K = N; K >= 1; K--) {
    if (K > M) {
      dtrmv('Lower', 'No transpose', 'Unit', M, AFAC, LDAFAC,
          AFAC(1, K).asArray(), 1);
    } else {
      // Compute elements (K+1:M,K)

      final T = AFAC[K][K];
      if (K + 1 <= M) {
        dscal(M - K, T, AFAC(K + 1, K).asArray(), 1);
        dgemv('No transpose', M - K, K - 1, ONE, AFAC(K + 1, 1), LDAFAC,
            AFAC(1, K).asArray(), 1, ONE, AFAC(K + 1, K).asArray(), 1);
      }

      // Compute the (K,K) element

      AFAC[K][K] = T +
          ddot(K - 1, AFAC(K, 1).asArray(), LDAFAC, AFAC(1, K).asArray(), 1);

      // Compute elements (1:K-1,K)

      dtrmv('Lower', 'No transpose', 'Unit', K - 1, AFAC, LDAFAC,
          AFAC(1, K).asArray(), 1);
    }
  }
  dlaswp(N, AFAC, LDAFAC, 1, min(M, N), IPIV, -1);

  // Compute the difference  L*U - A  and store in AFAC.

  for (var J = 1; J <= N; J++) {
    for (var I = 1; I <= M; I++) {
      AFAC[I][J] = AFAC[I][J] - A[I][J];
    }
  }

  // Compute norm( L*U - A ) / ( N * norm(A) * EPS )

  RESID.value = dlange('1', M, N, AFAC, LDAFAC, RWORK);

  if (ANORM <= ZERO) {
    if (RESID.value != ZERO) RESID.value = ONE / EPS;
  } else {
    RESID.value = ((RESID.value / N.toDouble()) / ANORM) / EPS;
  }
}
