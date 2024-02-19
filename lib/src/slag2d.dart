import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';

void slag2d(
  final int M,
  final int N,
  final Matrix<double> SA_,
  final int LDSA,
  final Matrix<double> A_,
  final int LDA,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final SA = SA_.dim(LDSA);
  final A = A_.dim(LDA);
  int I, J;

  INFO.value = 0;
  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      A[I][J] = SA[I][J];
    }
  }
}
