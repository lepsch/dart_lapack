import 'package:lapack/src/box.dart';

void dlabad(Box<double> SMALL, Box<double> LARGE) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // If it looks like we're on a Cray, take the square root of
  // SMALL and LARGE to avoid overflow and underflow problems.

  // if (log10(LARGE.value) > 2000.0) {
  //   SMALL.value = sqrt(SMALL.value);
  //   LARGE.value = sqrt(LARGE.value);
  // }
}
