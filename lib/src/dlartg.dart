import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/la_constants.dart';

const _zero = dzero;
const _one = done;
final _safmin = dsafmin;
final _safmax = dsafmax;

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
  final rtmin = sqrt(_safmin);
  final rtmax = sqrt(_safmax / 2);

  final f1 = f.abs();
  final g1 = g.abs();
  if (g == _zero) {
    c.value = _one;
    s.value = _zero;
    r.value = f;
  } else if (f == _zero) {
    c.value = _zero;
    s.value = sign(_one, g).toDouble();
    r.value = g1;
  } else if (f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax) {
    final d = sqrt(f * f + g * g);
    c.value = f1 / d;
    r.value = sign(d, f).toDouble();
    s.value = g / r.value;
  } else {
    final u = min(_safmax, max(_safmin, max(f1, g1)));
    final fs = f / u;
    final gs = g / u;
    final d = sqrt(fs * fs + gs * gs);
    c.value = fs.abs() / d;
    r.value = sign(d, f).toDouble();
    s.value = gs / r.value;
    r.value *= u;
  }
}
