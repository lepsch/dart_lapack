import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zlapmt(
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
  int I, II, IN, J;
  Complex TEMP;

  if (N <= 1) return;

  for (I = 1; I <= N; I++) {
    K[I] = -K[I];
  }

  if (FORWRD) {
    // Forward permutation

    forwardPermutation:
    for (I = 1; I <= N; I++) {
      if (K[I] > 0) continue forwardPermutation;

      J = I;
      K[J] = -K[J];
      IN = K[J];

      while (true) {
        if (K[IN] > 0) continue forwardPermutation;

        for (II = 1; II <= M; II++) {
          TEMP = X[II][J];
          X[II][J] = X[II][IN];
          X[II][IN] = TEMP;
        }

        K[IN] = -K[IN];
        J = IN;
        IN = K[IN];
      }
    }
  } else {
    // Backward permutation
    backwardPermutation:
    for (I = 1; I <= N; I++) {
      if (K[I] > 0) continue backwardPermutation;

      K[I] = -K[I];
      J = K[I];
      while (true) {
        if (J == I) continue backwardPermutation;

        for (II = 1; II <= M; II++) {
          TEMP = X[II][I];
          X[II][I] = X[II][J];
          X[II][J] = TEMP;
        }

        K[J] = -K[J];
        J = K[J];
      }
    }
  }
}
