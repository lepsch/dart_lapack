import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dlaqsp(
  final String UPLO,
  final int N,
  final Array<double> AP_,
  final Array<double> S_,
  final double SCOND,
  final double AMAX,
  final Box<String> EQUED,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.having();
  final S = S_.having();
  const ONE = 1.0, THRESH = 0.1;
  int I, J, JC;
  double CJ, LARGE, SMALL;

  // Quick return if possible

  if (N <= 0) {
    EQUED.value = 'N';
    return;
  }

  // Initialize LARGE and SMALL.

  SMALL = dlamch('Safe minimum') / dlamch('Precision');
  LARGE = ONE / SMALL;

  if (SCOND >= THRESH && AMAX >= SMALL && AMAX <= LARGE) {
    // No equilibration

    EQUED.value = 'N';
  } else {
    // Replace A by diag(S) * A * diag(S).

    if (lsame(UPLO, 'U')) {
      // Upper triangle of A is stored.

      JC = 1;
      for (J = 1; J <= N; J++) {
        CJ = S[J];
        for (I = 1; I <= J; I++) {
          AP[JC + I - 1] = CJ * S[I] * AP[JC + I - 1];
        }
        JC += J;
      }
    } else {
      // Lower triangle of A is stored.

      JC = 1;
      for (J = 1; J <= N; J++) {
        CJ = S[J];
        for (I = J; I <= N; I++) {
          AP[JC + I - J] = CJ * S[I] * AP[JC + I - J];
        }
        JC += N - J + 1;
      }
    }
    EQUED.value = 'Y';
  }
}
