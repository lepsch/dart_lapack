import 'dart:math';

const _log10e = 0.43429448190325182765;

double log10(final double x) {
  return (_log10e * log(x));
}
