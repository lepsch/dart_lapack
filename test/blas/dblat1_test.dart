import 'package:lapack/blas.dart';

import '../test_driver.dart';
import 'dblat1.dart';

void main() {
  Nout nout = NullNout();
  dblat1(nout, asyncTestDriver);
}
