N radix<N extends num>(final N n) {
  return switch (n) {
    double() => 2.0,
    int() => 2,
  } as N;
}
