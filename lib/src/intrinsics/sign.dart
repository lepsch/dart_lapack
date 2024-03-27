T sign<T extends num>(final T a, final num b) {
  return a.abs() * (b.isNegative ? -1 : 1) as T;
}
