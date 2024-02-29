import 'dart:math';

import 'package:lapack/src/intrinsics/digits.dart';
import 'package:lapack/src/intrinsics/radix.dart';

const double Function(double _) epsilon = epsilon64;

double epsilon32(final double _) => pow(radix32(_), 1 - digits32(_)).toDouble();
double epsilon64(final double _) => pow(radix(_), 1 - digits(_)).toDouble();
