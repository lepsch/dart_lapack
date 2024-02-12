import 'package:lapack/src/blas/lsame.dart';

bool lsamen(final int N, final String CA, final String CB) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  if (CA.length < N || CB.length < N) return false;

  // Do for each character in the two strings.

  for (var I = 1; I <= N; I++) {
    // Test if the characters are equal using lsame.
    if (!lsame(CA[I - 1], CB[I - 1])) return false;
  }
  return true;
}
