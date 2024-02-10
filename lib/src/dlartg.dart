import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/la_constants.dart';

const wp = dp;
const zero = dzero;
const one = done;
const half = dhalf;

void dlartg(
  final double f,
  final double g,
  final Box<double> c,
  final Box<double> s,
  final Box<double> r,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
//    February 2021

  double d, f1, fs, g1, gs, u, rtmin, rtmax;
  rtmin = sqrt(safmin);
  rtmax = sqrt(safmax / 2);

  f1 = f.abs();
  g1 = g.abs();
  if (g == zero) {
    c.value = one;
    s.value = zero;
    r.value = f;
  } else if (f == zero) {
    c.value = zero;
    s.value = sign(one, g).toDouble();
    r.value = g1;
  } else if (f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax) {
    d = sqrt(f * f + g * g);
    c.value = f1 / d;
    r.value = sign(d, f).toDouble();
    s.value = g / r.value;
  } else {
    u = min(safmax, max(safmin, max(f1, g1)));
    fs = f / u;
    gs = g / u;
    d = sqrt(fs * fs + gs * gs);
    c.value = fs.abs() / d;
    r.value = sign(d, f).toDouble();
    s.value = gs / r.value;
    r.value = r.value * u;
  }
}
