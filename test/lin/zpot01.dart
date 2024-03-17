import 'package:lapack/src/blas/zdotc.dart';
import 'package:lapack/src/blas/zher.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/blas/ztrmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlanhe.dart';

void zpot01(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> AFAC_,
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

  // Exit with RESID.value = 1/EPS if ANORM = 0.

  final EPS = dlamch('Epsilon');
  final ANORM = zlanhe('1', UPLO, N, A, LDA, RWORK);
  if (ANORM <= ZERO) {
    RESID.value = ONE / EPS;
    return;
  }

  // Check the imaginary parts of the diagonal elements and return with
  // an error code if any are nonzero.

  for (var J = 1; J <= N; J++) {
    if (AFAC[J][J].imaginary != ZERO) {
      RESID.value = ONE / EPS;
      return;
    }
  }

  // Compute the product U**H * U, overwriting U.

  if (lsame(UPLO, 'U')) {
    for (var K = N; K >= 1; K--) {
      // Compute the (K,K) element of the result.

      final TR =
          zdotc(K, AFAC(1, K).asArray(), 1, AFAC(1, K).asArray(), 1).real;
      AFAC[K][K] = TR.toComplex();

      // Compute the rest of column K.

      ztrmv('Upper', 'Conjugate', 'Non-unit', K - 1, AFAC, LDAFAC,
          AFAC(1, K).asArray(), 1);
    }

    // Compute the product L * L**H, overwriting L.
  } else {
    for (var K = N; K >= 1; K--) {
      // Add a multiple of column K of the factor L to each of
      // columns K+1 through N.

      if (K + 1 <= N) {
        zher('Lower', N - K, ONE, AFAC(K + 1, K).asArray(), 1,
            AFAC(K + 1, K + 1), LDAFAC);
      }

      // Scale column K by the diagonal element.

      final TC = AFAC[K][K];
      zscal(N - K + 1, TC, AFAC(K, K).asArray(), 1);
    }
  }

  // Compute the difference L * L**H - A (or U**H * U - A).

  if (lsame(UPLO, 'U')) {
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= J - 1; I++) {
        AFAC[I][J] -= A[I][J];
      }
      AFAC[J][J] -= A[J][J].real.toComplex();
    }
  } else {
    for (var J = 1; J <= N; J++) {
      AFAC[J][J] -= A[J][J].real.toComplex();
      for (var I = J + 1; I <= N; I++) {
        AFAC[I][J] -= A[I][J];
      }
    }
  }

  // Compute norm(L*U - A) / ( N * norm(A) * EPS )

  RESID.value = zlanhe('1', UPLO, N, AFAC, LDAFAC, RWORK);

  RESID.value = ((RESID.value / N) / ANORM) / EPS;
}
