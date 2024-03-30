import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dsyr.dart';
import 'package:lapack/src/blas/dtrmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';

void dpot01(
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> AFAC_,
  final int LDAFAC,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final AFAC = AFAC_.having(ld: LDAFAC);
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if N = 0.

  if (N <= 0) {
    RESID.value = ZERO;
    return;
  }

  // Exit with RESID = 1/EPS if ANORM = 0.

  final EPS = dlamch('Epsilon');
  final ANORM = dlansy('1', UPLO, N, A, LDA, RWORK);
  if (ANORM <= ZERO) {
    RESID.value = ONE / EPS;
    return;
  }

  // Compute the product U**T * U, overwriting U.

  if (lsame(UPLO, 'U')) {
    for (var K = N; K >= 1; K--) {
      // Compute the (K,K) element of the result.

      final T = ddot(K, AFAC(1, K).asArray(), 1, AFAC(1, K).asArray(), 1);
      AFAC[K][K] = T;

      // Compute the rest of column K.

      dtrmv('Upper', 'Transpose', 'Non-unit', K - 1, AFAC, LDAFAC,
          AFAC(1, K).asArray(), 1);
    }

    // Compute the product L * L**T, overwriting L.
  } else {
    for (var K = N; K >= 1; K--) {
      // Add a multiple of column K of the factor L to each of
      // columns K+1 through N.

      if (K + 1 <= N) {
        dsyr('Lower', N - K, ONE, AFAC(K + 1, K).asArray(), 1,
            AFAC(K + 1, K + 1), LDAFAC);
      }

      // Scale column K by the diagonal element.

      final T = AFAC[K][K];
      dscal(N - K + 1, T, AFAC(K, K).asArray(), 1);
    }
  }

  // Compute the difference L * L**T - A (or U**T * U - A).

  if (lsame(UPLO, 'U')) {
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= J; I++) {
        AFAC[I][J] -= A[I][J];
      }
    }
  } else {
    for (var J = 1; J <= N; J++) {
      for (var I = J; I <= N; I++) {
        AFAC[I][J] -= A[I][J];
      }
    }
  }

  // Compute norm(L*U - A) / ( N * norm(A) * EPS )

  RESID.value = dlansy('1', UPLO, N, AFAC, LDAFAC, RWORK);

  RESID.value = ((RESID.value / N) / ANORM) / EPS;
}
