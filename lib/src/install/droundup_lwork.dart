import 'package:lapack/src/intrinsics/epsilon.dart';

double droundup_lwork(int LWORK) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  var result = LWORK.toDouble();

  if (result.toInt() < LWORK) {
    // Force round up of LWORK
    result *= (1.0 + epsilon(0.0));
  }
  return result;
}
