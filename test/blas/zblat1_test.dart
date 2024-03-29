import 'package:lapack/src/nio.dart';

import '../test_driver.dart';
import 'zblat1.dart';

void main() {
  Nout nout = NullNout();
  zblat1(nout, asyncTestDriver);
}
