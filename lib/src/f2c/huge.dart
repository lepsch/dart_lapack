N huge<N extends num>(final N n) {
  return switch (n) {
    double() => double.maxFinite,
    int() => 0x7fffffffffffffff,
  } as N;
}
