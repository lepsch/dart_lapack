import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/ilaver.dart';

void main() {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  final MAJOR = Box(0), MINOR = Box(0), PATCH = Box(0);

  ilaver(MAJOR, MINOR, PATCH);
  print('LAPACK $MAJOR.$MINOR.$PATCH');
}
