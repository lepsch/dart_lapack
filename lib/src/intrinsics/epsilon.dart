import 'dart:math';

import 'package:lapack/src/intrinsics/digits.dart';
import 'package:lapack/src/intrinsics/radix.dart';

double epsilon(final double _) => pow(radix(_), 1 - digits(_)).toDouble();
