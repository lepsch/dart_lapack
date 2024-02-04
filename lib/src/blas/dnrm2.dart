import 'dart:math';

import 'package:lapack/src/matrix.dart';

double dnrm2(final int n, Array<double> x, final int incx) {
// -- Reference BLAS level1 routine (version 3.9.1) --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
//    March 2021

  // .. Constants ..
  const zero = 0.0;
  const one = 1.0;
  const maxN = huge(0.0);

  // .. Blue's scaling constants ..
  //  const tsml = pow(real(radix(0.0), wp),ceiling(
  //      (minexponent(0.0) - 1) * 0.5));
  //  const tbig = pow(real(radix(0.0), wp),floor(
  //      (maxexponent(0.0) - digits(0.0) + 1) * 0.5));
  //  const ssml = pow(real(radix(0.0), wp),( - floor(
  //      (minexponent(0.0) - digits(0.0)) * 0.5)));
  //  const sbig = pow(real(radix(0.0), wp),( - ceiling(
  //      (maxexponent(0.0) + digits(0.0) - 1) * 0.5)));

  int i, ix;
  bool notbig;
  double abig, amed, asml, ax = 0, scl, sumsq, ymax = 0, ymin;

  // Quick return if possible

  if (n <= 0) return zero;

  scl = one;
  sumsq = zero;

  // Compute the sum of squares in 3 accumulators:
  //    abig -- sums of squares scaled down to avoid overflow
  //    asml -- sums of squares scaled up to avoid underflow
  //    amed -- sums of squares that do not require scaling
  // The thresholds and multipliers are
  //    tbig -- values bigger than this are scaled down by sbig
  //    tsml -- values smaller than this are scaled up by ssml

  notbig = true;
  asml = zero;
  amed = zero;
  abig = zero;
  ix = 1;
  if (incx < 0) ix = 1 - (n - 1) * incx;
  for (i = 1; i <= n; i++) {
    ax = x(ix).abs();
    if (ax > tbig) {
      abig = abig + pow((ax * sbig), 2);
      notbig = false;
    } else if (ax < tsml) {
      if (notbig) asml = asml + pow((ax * ssml), 2);
    } else {
      amed = amed + pow(ax, 2);
    }
    ix = ix + incx;
  }

  // Combine abig and amed or amed and asml if more than one
  // accumulator was used.

  if (abig > zero) {
    // Combine abig and amed if abig > 0.

    if ((amed > zero) || (amed > maxN) || amed.isNaN) {
      abig = abig + (amed * sbig) * sbig;
    }
    scl = one / sbig;
    sumsq = abig;
  } else if (asml > zero) {
    // Combine amed and asml if asml > 0.

    if ((amed > zero) || (amed > maxN) || amed.isNaN) {
      amed = sqrt(amed);
      asml = sqrt(asml) / ssml;
      if (asml > amed) {
        ymin = amed;
        ymax = asml;
      } else {
        ymin = asml;
        ymax = amed;
      }
      scl = one;
      sumsq = pow(ymax, 2) * (one + pow((ymin / ymax), 2));
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
