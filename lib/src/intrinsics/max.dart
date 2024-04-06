import 'dart:math' as math;

T max<T extends num>(T a, T b) {
  if (a.isNaN) return b;
  if (b.isNaN) return a;
  return math.max(a, b);
}
