import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void zlaqge(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> R_,
  final Array<double> C_,
  final double ROWCND,
  final double COLCND,
  final double AMAX,
  final Box<String> EQUED,
) {
  final A = A_.having(ld: LDA);
  final R = R_.having();
  final C = C_.having();

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0, THRESH = 0.1;
  int I, J;
  double CJ, LARGE, SMALL;

  // Quick return if possible

  if (M <= 0 || N <= 0) {
    EQUED.value = 'N';
    return;
  }

  // Initialize LARGE and SMALL.

  SMALL = dlamch('Safe minimum') / dlamch('Precision');
  LARGE = ONE / SMALL;

  if (ROWCND >= THRESH && AMAX >= SMALL && AMAX <= LARGE) {
    // No row scaling

    if (COLCND >= THRESH) {
      // No column scaling

      EQUED.value = 'N';
    } else {
      // Column scaling

      for (J = 1; J <= N; J++) {
        // 20
        CJ = C[J];
        for (I = 1; I <= M; I++) {
          // 10
          A[I][J] = CJ.toComplex() * A[I][J];
        } // 10
      } // 20
      EQUED.value = 'C';
    }
  } else if (COLCND >= THRESH) {
    // Row scaling, no column scaling

    for (J = 1; J <= N; J++) {
      // 40
      for (I = 1; I <= M; I++) {
        // 30
        A[I][J] = R[I].toComplex() * A[I][J];
      } // 30
    } // 40
    EQUED.value = 'R';
  } else {
    // Row and column scaling

    for (J = 1; J <= N; J++) {
      // 60
      CJ = C[J];
      for (I = 1; I <= M; I++) {
        // 50
        A[I][J] = (CJ * R[I]).toComplex() * A[I][J];
      } // 50
    } // 60
    EQUED.value = 'B';
  }
}
