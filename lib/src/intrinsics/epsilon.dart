import 'dart:math';

import 'package:lapack/src/intrinsics/digits.dart';
import 'package:lapack/src/intrinsics/radix.dart';

const double Function(double v) epsilon = epsilon64;

double epsilon32(final double v) => pow(radix32(v), 1 - digits32(v)).toDouble();
double epsilon64(final double v) => pow(radix64(v), 1 - digits64(v)).toDouble();
