import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/intrinsics/maxexponent.dart';
import 'package:lapack/src/intrinsics/minexponent.dart';
import 'package:lapack/src/intrinsics/radix.dart';
import 'package:lapack/src/intrinsics/sign.dart';

void drotg(
  final Box<double> a,
  final Box<double> b,
  final Box<double> c,
  final Box<double> s,
) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  const zero = 0.0, one = 1.0;

  final safmin =
      pow(radix(0.0), max(minexponent(0.0) - 1, 1 - maxexponent(0.0)));
  final safmax =
      pow(radix(0.0), max(1 - minexponent(0.0), maxexponent(0.0) - 1));

  final anorm = a.value.abs();
  final bnorm = b.value.abs();
  if (bnorm == zero) {
    c.value = one;
    s.value = zero;
    b.value = zero;
  } else if (anorm == zero) {
    c.value = zero;
    s.value = one;
    a.value = b.value;
    b.value = one;
  } else {
    final scl = min(safmax, max(safmin, max(anorm, bnorm)));
    final sigma = sign(one, anorm > bnorm ? a.value : b.value);

    final r =
        sigma * (scl * sqrt(pow(a.value / scl, 2) + pow(b.value / scl, 2)));
    c.value = a.value / r;
    s.value = b.value / r;
    final double z;
    if (anorm > bnorm) {
      z = s.value;
    } else if (c.value != zero) {
      z = one / c.value;
    } else {
      z = one;
    }
    a.value = r;
    b.value = z;
  }
}
