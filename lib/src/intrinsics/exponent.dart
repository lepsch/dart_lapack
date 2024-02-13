import 'dart:typed_data';

import 'package:lapack/src/intrinsics/maxexponent.dart';

int exponent(final double x) {
  final bias = maxexponent(x) - 1;
  return ((Float64List.fromList([x]).buffer.asUint32List(0).last &
              0x7fffffff) >>>
          (52 - 32)) -
      bias + 1;
}
