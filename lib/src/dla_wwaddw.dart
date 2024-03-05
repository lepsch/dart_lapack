import 'package:lapack/src/matrix.dart';

void dla_wwaddw(
  final int N,
  final Array<double> X_,
  final Array<double> Y_,
  final Array<double> W_,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having();
  final Y = Y_.having();
  final W = W_.having();
  double S;

  for (var I = 1; I <= N; I++) {
    S = X[I] + W[I];
    S = (S + S) - S;
    Y[I] = ((X[I] - S) + W[I]) + Y[I];
    X[I] = S;
  }
}
