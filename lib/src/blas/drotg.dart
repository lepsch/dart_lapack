import 'dart:math';

import 'package:lapack/src/box.dart';
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

  const zero = 0.0;
  const one = 1.0;
  // const wp = 1.0;
  //  double safmin = real(radix(0._wp),wp)**max(
  //     minexponent(0._wp)-1,
  //     1-maxexponent(0._wp)
  //  );
  //  double safmax = real(radix(0._wp),wp)**max(
  //     1-minexponent(0._wp),
  //     maxexponent(0._wp)-1
  //  );
  double anorm, bnorm, scl, sigma, r, z;

  anorm = a.value.abs();
  bnorm = b.value.abs();
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
    scl = min(safmax, max(safmin, max(anorm, bnorm)));
    if (anorm > bnorm) {
      sigma = sign(one, a.value).toDouble();
    } else {
      sigma = sign(one, b.value).toDouble();
    }
    r = sigma * (scl * sqrt(pow((a.value / scl), 2) + pow((b.value / scl), 2)));
    c.value = a.value / r;
    s.value = b.value / r;
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
