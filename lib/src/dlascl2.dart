import 'package:lapack/src/matrix.dart';

void dlascl2(
  final int M,
  final int N,
  final Array<double> D,
  final Matrix<double> X,
  final int LDX,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  int I, J;

  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      X[I][J] = X[I][J] * D[I];
    }
  }
}
