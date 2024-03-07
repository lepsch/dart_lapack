import 'dart:math';

import 'package:lapack/src/matrix.dart';

bool dgennd(final int M, final int N, final Matrix<double> A_, final int LDA) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  const ZERO = 0.0;

  final K = min(M, N);
  for (var I = 1; I <= K; I++) {
    if (A[I][I] < ZERO) {
      return false;
    }
  }
  return true;
}
