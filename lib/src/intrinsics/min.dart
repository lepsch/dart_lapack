import 'dart:math' as math;

T min<T extends num>(T a, T b) {
  if (a.isNaN) return b;
  if (b.isNaN) return a;
  return math.min(a, b);
}
