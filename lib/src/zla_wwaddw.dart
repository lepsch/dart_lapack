import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zla_wwaddw(
  final int N,
  final Array<Complex> X_,
  final Array<Complex> Y_,
  final Array<Complex> W_,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having();
  final Y = Y_.having();
  final W = W_.having();

  for (var I = 1; I <= N; I++) {
    var S = X[I] + W[I];
    S = (S + S) - S;
    Y[I] = ((X[I] - S) + W[I]) + Y[I];
    X[I] = S;
  }
}
