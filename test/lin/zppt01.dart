import 'package:lapack/src/blas/zdotc.dart';
import 'package:lapack/src/blas/zhpr.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/blas/ztpmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlanhp.dart';

void zppt01(
  final String UPLO,
  final int N,
  final Array<Complex> A_,
  final Array<Complex> AFAC_,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having();
  final AFAC = AFAC_.having();
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if N = 0

  if (N <= 0) {
    RESID.value = ZERO;
    return;
  }

  // Exit with RESID = 1/EPS if ANORM = 0.

  final EPS = dlamch('Epsilon');
  final ANORM = zlanhp('1', UPLO, N, A, RWORK);
  if (ANORM <= ZERO) {
    RESID.value = ONE / EPS;
    return;
  }

  // Check the imaginary parts of the diagonal elements and return with
  // an error code if any are nonzero.

  var KC = 1;
  if (lsame(UPLO, 'U')) {
    for (var K = 1; K <= N; K++) {
      if (AFAC[KC].imaginary != ZERO) {
        RESID.value = ONE / EPS;
        return;
      }
      KC += K + 1;
    }
  } else {
    for (var K = 1; K <= N; K++) {
      if (AFAC[KC].imaginary != ZERO) {
        RESID.value = ONE / EPS;
        return;
      }
      KC += N - K + 1;
    }
  }

  // Compute the product U'*U, overwriting U.

  if (lsame(UPLO, 'U')) {
    KC = (N * (N - 1)) ~/ 2 + 1;
    for (var K = N; K >= 1; K--) {
      // Compute the (K,K) element of the result.

      final TR = zdotc(K, AFAC(KC), 1, AFAC(KC), 1).real;
      AFAC[KC + K - 1] = TR.toComplex();

      // Compute the rest of column K.

      if (K > 1) {
        ztpmv('Upper', 'Conjugate', 'Non-unit', K - 1, AFAC, AFAC(KC), 1);
        KC -= (K - 1);
      }
    }

    // Compute the difference  L*L' - A

    KC = 1;
    for (var K = 1; K <= N; K++) {
      for (var I = 1; I <= K - 1; I++) {
        AFAC[KC + I - 1] -= A[KC + I - 1];
      }
      AFAC[KC + K - 1] -= A[KC + K - 1].real.toComplex();
      KC += K;
    }

    // Compute the product L*L', overwriting L.
  } else {
    KC = (N * (N + 1)) ~/ 2;
    for (var K = N; K >= 1; K--) {
      // Add a multiple of column K of the factor L to each of
      // columns K+1 through N.

      if (K < N) {
        zhpr('Lower', N - K, ONE, AFAC(KC + 1), 1, AFAC(KC + N - K + 1));
      }

      // Scale column K by the diagonal element.

      final TC = AFAC[KC];
      zscal(N - K + 1, TC, AFAC(KC), 1);

      KC -= (N - K + 2);
    }

    // Compute the difference  U'*U - A

    KC = 1;
    for (var K = 1; K <= N; K++) {
      AFAC[KC] -= A[KC].real.toComplex();
      for (var I = K + 1; I <= N; I++) {
        AFAC[KC + I - K] -= A[KC + I - K];
      }
      KC += N - K + 1;
    }
  }

  // Compute norm( L*U - A ) / ( N * norm(A) * EPS )

  RESID.value = zlanhp('1', UPLO, N, AFAC, RWORK);

  RESID.value = ((RESID.value / N) / ANORM) / EPS;
}
