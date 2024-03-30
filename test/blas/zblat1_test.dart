import 'package:lapack/blas.dart';

import '../test_driver.dart';
import 'zblat1.dart';

void main() {
  Nout nout = NullNout();
  zblat1(nout, asyncTestDriver);
}
