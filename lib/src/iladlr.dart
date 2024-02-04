import 'dart:math';

import 'package:lapack/src/matrix.dart';

int iladlr(final int M, final int N, final Matrix<double> A, final int LDA) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;
  int I, J;

  // Quick test for the common case where one corner is non-zero.
  if (M == 0) {
    return M;
  }

  if (A(M, 1) != ZERO || A(M, N) != ZERO) {
    return M;
  }

  // Scan up each column tracking the last zero row seen.
  var result = 0;
  for (J = 1; J <= N; J++) {
    I = M;
    while ((A[max(I, 1)][J] == ZERO) && (I >= 1)) {
      I = I - 1;
    }
    result = max(result, I);
  }
  return result;
}
