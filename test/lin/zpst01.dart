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

void zpst01(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> AFAC_,
  final int LDAFAC,
  final Matrix<Complex> PERM_,
  final int LDPERM,
  final Array<int> PIV_,
  final Array<double> RWORK_,
  final Box<double> RESID,
  final int RANK,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final PIV = PIV_.having();
  final A = A_.having(ld: LDA);
  final AFAC = AFAC_.having(ld: LDAFAC);
  final PERM = PERM_.having(ld: LDPERM);
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if N = 0.

  if (N <= 0) {
    RESID.value = ZERO;
    return;
  }

  // Exit with RESID = 1/EPS if ANORM = 0.

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

  // Compute the product U'*U, overwriting U.

  if (lsame(UPLO, 'U')) {
    if (RANK < N) {
      for (var J = RANK + 1; J <= N; J++) {
        for (var I = RANK + 1; I <= J; I++) {
          AFAC[I][J] = Complex.zero;
        }
      }
    }

    for (var K = N; K >= 1; K--) {
      // Compute the (K,K) element of the result.

      final TR =
          zdotc(K, AFAC(1, K).asArray(), 1, AFAC(1, K).asArray(), 1).real;
      AFAC[K][K] = TR.toComplex();

      // Compute the rest of column K.

      ztrmv('Upper', 'Conjugate', 'Non-unit', K - 1, AFAC, LDAFAC,
          AFAC(1, K).asArray(), 1);
    }

    // Compute the product L*L', overwriting L.
  } else {
    if (RANK < N) {
      for (var J = RANK + 1; J <= N; J++) {
        for (var I = J; I <= N; I++) {
          AFAC[I][J] = Complex.zero;
        }
      }
    }

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

  // Form P*L*L'*P' or P*U'*U*P'

  if (lsame(UPLO, 'U')) {
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= N; I++) {
        if (PIV[I] <= PIV[J]) {
          if (I <= J) {
            PERM[PIV[I]][PIV[J]] = AFAC[I][J];
          } else {
            PERM[PIV[I]][PIV[J]] = AFAC[J][I].conjugate();
          }
        }
      }
    }
  } else {
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= N; I++) {
        if (PIV[I] >= PIV[J]) {
          if (I >= J) {
            PERM[PIV[I]][PIV[J]] = AFAC[I][J];
          } else {
            PERM[PIV[I]][PIV[J]] = AFAC[J][I].conjugate();
          }
        }
      }
    }
  }

  // Compute the difference  P*L*L'*P' - A (or P*U'*U*P' - A).

  if (lsame(UPLO, 'U')) {
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= J - 1; I++) {
        PERM[I][J] -= A[I][J];
      }
      PERM[J][J] -= A[J][J].real.toComplex();
    }
  } else {
    for (var J = 1; J <= N; J++) {
      PERM[J][J] -= A[J][J].real.toComplex();
      for (var I = J + 1; I <= N; I++) {
        PERM[I][J] -= A[I][J];
      }
    }
  }

  // Compute norm( P*L*L'P - A ) / ( N * norm(A) * EPS ), or
  // ( P*U'*U*P' - A )/ ( N * norm(A) * EPS ).

  RESID.value = zlanhe('1', UPLO, N, PERM, LDAFAC, RWORK);

  RESID.value = ((RESID.value / N) / ANORM) / EPS;
}
