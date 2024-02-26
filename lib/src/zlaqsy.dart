import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void zlaqsy(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> S_,
  final double SCOND,
  final double AMAX,
  final Box<String> EQUED,
) {
  final A = A_.dim(LDA);
  final S = S_.dim();

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0, THRESH = 0.1;
  int I, J;
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

      for (J = 1; J <= N; J++) {
        // 20
        CJ = S[J];
        for (I = 1; I <= J; I++) {
          // 10
          A[I][J] = (CJ * S[I]).toComplex() * A[I][J];
        } // 10
      } // 20
    } else {
      // Lower triangle of A is stored.

      for (J = 1; J <= N; J++) {
        // 40
        CJ = S[J];
        for (I = J; I <= N; I++) {
          // 30
          A[I][J] = (CJ * S[I]).toComplex() * A[I][J];
        } // 30
      } // 40
    }
    EQUED.value = 'Y';
  }
}
