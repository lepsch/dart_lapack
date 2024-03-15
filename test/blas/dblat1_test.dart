import 'package:lapack/src/nio.dart';

import '../test_driver.dart';
import 'dblat1.dart';

void main() {
  Nout nout = NullNout();
  dblat1(nout, dartTestDriver);
}
