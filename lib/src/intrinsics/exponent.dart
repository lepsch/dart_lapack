import 'dart:typed_data';

import 'package:lapack/src/intrinsics/digits.dart';
import 'package:lapack/src/intrinsics/maxexponent.dart';

const int Function(double x) exponent = exponent64;

int exponent64(final double x) {
  final bias = maxexponent(x) - 1;
  final d = digits(x) - 1;
  return ((Float64List.fromList([x]).buffer.asUint32List(0).last &
              0x7fffffff) >>>
          (d - 32)) -
      bias +
      1;
}

int exponent32(final double x) {
  final bias = maxexponent32(x) - 1;
  final d = digits32(x) - 1;
  return ((Float64List.fromList([x]).buffer.asUint32List(0).last &
              0x7fffffff) >>>
          (d - 32)) -
      bias +
      1;
}
