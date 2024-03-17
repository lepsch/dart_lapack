import 'dart:math';

import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dsyr.dart';
import 'package:lapack/src/blas/dtrmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlansb.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';

void dpbt01(
  final String UPLO,
  final int N,
  final int KD,
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
  final ANORM = dlansb('1', UPLO, N, KD, A, LDA, RWORK);
  if (ANORM <= ZERO) {
    RESID.value = ONE / EPS;
    return;
  }

  // Compute the product U'*U, overwriting U.

  if (lsame(UPLO, 'U')) {
    for (var K = N; K >= 1; K--) {
      final KC = max(1, KD + 2 - K);
      final KLEN = KD + 1 - KC;

      // Compute the (K,K) element of the result.

      final T =
          ddot(KLEN + 1, AFAC(KC, K).asArray(), 1, AFAC(KC, K).asArray(), 1);
      AFAC[KD + 1][K] = T;

      // Compute the rest of column K.

      if (KLEN > 0)
        dtrmv('Upper', 'Transpose', 'Non-unit', KLEN, AFAC(KD + 1, K - KLEN),
            LDAFAC - 1, AFAC(KC, K).asArray(), 1);
    }

    // UPLO = 'L':  Compute the product L*L', overwriting L.
  } else {
    for (var K = N; K >= 1; K--) {
      final KLEN = min(KD, N - K);

      // Add a multiple of column K of the factor L to each of
      // columns K+1 through N.

      if (KLEN > 0) {
        dsyr('Lower', KLEN, ONE, AFAC(2, K).asArray(), 1, AFAC(1, K + 1),
            LDAFAC - 1);
      }

      // Scale column K by the diagonal element.

      final T = AFAC[1][K];
      dscal(KLEN + 1, T, AFAC(1, K).asArray(), 1);
    }
  }

  // Compute the difference  L*L' - A  or  U'*U - A.

  if (lsame(UPLO, 'U')) {
    for (var J = 1; J <= N; J++) {
      final MU = max(1, KD + 2 - J);
      for (var I = MU; I <= KD + 1; I++) {
        AFAC[I][J] -= A[I][J];
      }
    }
  } else {
    for (var J = 1; J <= N; J++) {
      final ML = min(KD + 1, N - J + 1);
      for (var I = 1; I <= ML; I++) {
        AFAC[I][J] -= A[I][J];
      }
    }
  }

  // Compute norm( L*L' - A ) / ( N * norm(A) * EPS )

  RESID.value = dlansb('I', UPLO, N, KD, AFAC, LDAFAC, RWORK);

  RESID.value = ((RESID.value / N) / ANORM) / EPS;
}
