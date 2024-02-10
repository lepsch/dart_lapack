import 'dart:math';

import 'package:lapack/src/f2c/digits.dart';
import 'package:lapack/src/f2c/radix.dart';

double epsilon(final double _) => pow(radix(_), 1 - digits(_)).toDouble();
