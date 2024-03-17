import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zlarscl2(
  final int M,
  final int N,
  final Array<double> D_,
  final Matrix<Complex> X_,
  final int LDX,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final X = X_.having(ld: LDX);
  int I, J;

  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      X[I][J] /= D[I].toComplex();
    }
  }
}
