import 'dart:typed_data';

int exponent(final double x) {
  return Float64List.fromList([x]).buffer.asUint32List(0).first >> (23 - 16);
}
