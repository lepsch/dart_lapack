import 'dart:math';

import 'package:lapack/src/box.dart';

final _rnd = Random();

void random_number(final Box<double> r) {
  r.value = _rnd.nextDouble();
}
