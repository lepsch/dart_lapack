import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/slamch.dart';
import 'package:lapack/src/matrix.dart';

void dlag2s(
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> SA_,
  final int LDSA,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final SA = SA_.having(ld: LDSA);

  final RMAX = slamch('O');
  for (var J = 1; J <= N; J++) {
    for (var I = 1; I <= M; I++) {
      if ((A[I][J] < -RMAX) || (A[I][J] > RMAX)) {
        INFO.value = 1;
        return;
      }
      SA[I][J] = A[I][J];
    }
  }
  INFO.value = 0;
}
