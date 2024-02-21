import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zlat2c(
  final String UPLO,
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
  final A = A_.dim(LDA);
  final SA = SA_.dim(LDSA);
  int I, J;
  double RMAX;
  bool UPPER;

  RMAX = slamch('O');
  UPPER = lsame(UPLO, 'U');
  if (UPPER) {
    for (J = 1; J <= N; J++) {
      // 20
      for (I = 1; I <= J; I++) {
        // 10
        if (((A[I][J]).toDouble() < -RMAX) ||
            ((A[I][J]).toDouble() > RMAX) ||
            (A[I][J].imaginary < -RMAX) ||
            (A[I][J].imaginary > RMAX)) {
          INFO.value = 1;
          return;
        }
        SA[I][J] = A[I][J];
      } // 10
    } // 20
  } else {
    for (J = 1; J <= N; J++) {
      // 40
      for (I = J; I <= N; I++) {
        // 30
        if (((A[I][J]).toDouble() < -RMAX) ||
            ((A[I][J]).toDouble() > RMAX) ||
            (A[I][J].imaginary < -RMAX) ||
            (A[I][J].imaginary > RMAX)) {
          INFO.value = 1;
          return;
        }
        SA[I][J] = A[I][J];
      } // 30
    } // 40
  }
}
