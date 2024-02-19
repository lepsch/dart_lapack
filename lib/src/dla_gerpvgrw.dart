import 'dart:math';

import 'package:lapack/src/matrix.dart';

double dla_gerpvgrw(
  final int N,
  final int NCOLS,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> AF_,
  final int LDAF,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final AF = AF_.dim(LDAF);
  int I, J;
  double AMAX, UMAX, RPVGRW;

  RPVGRW = 1.0;

  for (J = 1; J <= NCOLS; J++) {
    AMAX = 0.0;
    UMAX = 0.0;
    for (I = 1; I <= N; I++) {
      AMAX = max(A[I][J].abs(), AMAX);
    }
    for (I = 1; I <= J; I++) {
      UMAX = max(AF[I][J].abs(), UMAX);
    }
    if (UMAX != 0.0) {
      RPVGRW = min(AMAX / UMAX, RPVGRW);
    }
  }
  return RPVGRW;
}
