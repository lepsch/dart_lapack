import 'dart:math';

import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

double zla_gerpvgrw(
  final int N,
  final int NCOLS,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> AF_,
  final int LDAF,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDAF);
  int I, J;
  double AMAX, UMAX, RPVGRW;

  RPVGRW = 1.0;

  for (J = 1; J <= NCOLS; J++) {
    AMAX = 0.0;
    UMAX = 0.0;
    for (I = 1; I <= N; I++) {
      AMAX = max(A[I][J].cabs1(), AMAX);
    }
    for (I = 1; I <= J; I++) {
      UMAX = max(AF[I][J].cabs1(), UMAX);
    }
    if (UMAX != 0.0) {
      RPVGRW = min(AMAX / UMAX, RPVGRW);
    }
  }
  return RPVGRW;
}
