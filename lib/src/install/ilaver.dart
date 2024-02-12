import 'package:lapack/src/box.dart';

void ilaver(final Box<int> VERS_MAJOR, final Box<int> VERS_MINOR, final Box<int> VERS_PATCH) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  VERS_MAJOR.value = 3;
  VERS_MINOR.value = 12;
  VERS_PATCH.value = 0;
}
