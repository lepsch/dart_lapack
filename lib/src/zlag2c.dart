import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/slamch.dart';
import 'package:lapack/src/matrix.dart';

void zlag2c(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> SA_,
  final int LDSA,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final SA = SA_.having(ld: LDSA);
  int I, J;
  double RMAX;

  RMAX = slamch('O');
  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      if ((A[I][J].real < -RMAX) ||
          (A[I][J].real > RMAX) ||
          (A[I][J].imaginary < -RMAX) ||
          (A[I][J].imaginary > RMAX)) {
        INFO.value = 1;
        return;
      }
      SA[I][J] = A[I][J];
    }
  }
  INFO.value = 0;
}
