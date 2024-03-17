import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void clag2z(
  final int M,
  final int N,
  final Matrix<Complex> SA_,
  final int LDSA,
  final Matrix<Complex> A_,
  final int LDA,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final SA = SA_.having(ld: LDSA);
  final A = A_.having(ld: LDA);
  int I, J;

  INFO.value = 0;
  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      A[I][J] = SA[I][J];
    }
  }
}
