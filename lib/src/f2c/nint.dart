int nint(final double x) {
  return (x >= 0 ? (x + .5).floor() : -(.5 - x).floor()).toInt();
}
