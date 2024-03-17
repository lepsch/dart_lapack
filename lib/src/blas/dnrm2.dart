import 'dart:math';

import 'package:lapack/src/intrinsics/digits.dart';
import 'package:lapack/src/intrinsics/huge.dart';
import 'package:lapack/src/intrinsics/maxexponent.dart';
import 'package:lapack/src/intrinsics/minexponent.dart';
import 'package:lapack/src/intrinsics/radix.dart';
import 'package:lapack/src/matrix.dart';

double dnrm2(final int n, final Array<double> x_, final int incx) {
// -- Reference BLAS level1 routine (version 3.9.1) --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final x = x_.having();
  const zero = 0.0;
  const one = 1.0;
  final maxN = huge(0.0);

  // .. Blue's scaling constants ..
  final tsml = pow(radix(0.0), ((minexponent(0.0) - 1) * 0.5).ceil());
  final tbig =
      pow(radix(0.0), ((maxexponent(0.0) - digits(0.0) + 1) * 0.5).floor());
  final ssml =
      pow(radix(0.0), -((minexponent(0.0) - digits(0.0)) * 0.5).floor());
  final sbig =
      pow(radix(0.0), -((maxexponent(0.0) + digits(0.0) - 1) * 0.5).ceil());

  // Quick return if possible

  if (n <= 0) return zero;

  // Compute the sum of squares in 3 accumulators:
  //    abig -- sums of squares scaled down to avoid overflow
  //    asml -- sums of squares scaled up to avoid underflow
  //    amed -- sums of squares that do not require scaling
  // The thresholds and multipliers are
  //    tbig -- values bigger than this are scaled down by sbig
  //    tsml -- values smaller than this are scaled up by ssml

  var notbig = true;
  var asml = zero;
  var amed = zero;
  var abig = zero;
  var ix = incx < 0 ? 1 - (n - 1) * incx : 1;
  for (var i = 1; i <= n; i++) {
    final ax = x[ix].abs();
    if (ax > tbig) {
      abig += pow(ax * sbig, 2);
      notbig = false;
    } else if (ax < tsml) {
      if (notbig) asml += pow(ax * ssml, 2);
    } else {
      amed += pow(ax, 2);
    }
    ix += incx;
  }

  // Combine abig and amed or amed and asml if more than one
  // accumulator was used.

  final double scl, sumsq;
  if (abig > zero) {
    // Combine abig and amed if abig > 0.

    if ((amed > zero) || (amed > maxN) || amed.isNaN) {
      abig += (amed * sbig) * sbig;
    }
    scl = one / sbig;
    sumsq = abig;
  } else if (asml > zero) {
    // Combine amed and asml if asml > 0.

    if ((amed > zero) || (amed > maxN) || amed.isNaN) {
      amed = sqrt(amed);
      asml = sqrt(asml) / ssml;
      final (ymin, ymax) = asml > amed ? (amed, asml) : (asml, amed);
      scl = one;
      sumsq = pow(ymax, 2) * (one + pow(ymin / ymax, 2));
    } else {
      scl = one / ssml;
      sumsq = asml;
    }
  } else {
    // Otherwise all values are mid-range

    scl = one;
    sumsq = amed;
  }
  return scl * sqrt(sumsq);
}
