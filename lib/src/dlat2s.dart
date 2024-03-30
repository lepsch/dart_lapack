import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/slamch.dart';
import 'package:lapack/src/matrix.dart';

void dlat2s(
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> SA_,
  final int LDSA,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final SA = SA_.having(ld: LDSA);

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  int I, J;
  double RMAX;
  bool UPPER;

  RMAX = slamch('O');
  UPPER = lsame(UPLO, 'U');
  if (UPPER) {
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= J; I++) {
        if ((A[I][J] < -RMAX) || (A[I][J] > RMAX)) {
          INFO.value = 1;
          return;
        }
        SA[I][J] = A[I][J];
      }
    }
  } else {
    for (J = 1; J <= N; J++) {
      for (I = J; I <= N; I++) {
        if ((A[I][J] < -RMAX) || (A[I][J] > RMAX)) {
          INFO.value = 1;
          return;
        }
        SA[I][J] = A[I][J];
      }
    }
  }
}
