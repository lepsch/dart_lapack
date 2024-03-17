import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zlapmr(
  final bool FORWRD,
  final int M,
  final int N,
  final Matrix<Complex> X_,
  final int LDX,
  final Array<int> K_,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having(ld: LDX);
  final K = K_.having();
  int I, IN = 0, J, JJ;
  Complex TEMP;

  if (M <= 1) return;

  for (I = 1; I <= M; I++) {
    K[I] = -K[I];
  }

  if (FORWRD) {
    // Forward permutation

    for (I = 1; I <= M; I++) {
      if (K[I] <= 0) {
        J = I;
        K[J] = -K[J];
        IN = K[J];
        while (K[IN] <= 0) {
          for (JJ = 1; JJ <= N; JJ++) {
            TEMP = X[J][JJ];
            X[J][JJ] = X[IN][JJ];
            X[IN][JJ] = TEMP;
          }

          K[IN] = -K[IN];
          J = IN;
          IN = K[IN];
        }
      }
    }
  } else {
    // Backward permutation

    for (I = 1; I <= M; I++) {
      if (K[I] <= 0) {
        K[I] = -K[I];
        J = K[I];
        while (J != I) {
          for (JJ = 1; JJ <= N; JJ++) {
            TEMP = X[I][JJ];
            X[I][JJ] = X[J][JJ];
            X[J][JJ] = TEMP;
          }

          K[J] = -K[J];
          J = K[J];
        }
      }
    }
  }
}
