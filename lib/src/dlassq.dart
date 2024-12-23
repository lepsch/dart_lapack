// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/disnan.dart';
import 'package:dart_lapack/src/la_constants.dart';
import 'package:dart_lapack/src/matrix.dart';

const wp = dp;
const zero = dzero;
const one = done;
final sbig = dsbig;
final ssml = dssml;
final tbig = dtbig;
final tsml = dtsml;

void dlassq(
  final int n,
  final Array<double> x_,
  final int incx,
  final Box<double> scale,
  final Box<double> sumsq,
) {
  final x = x_.having();
  //  use LA_XISNAN

  int i, ix;
  bool notbig;
  double abig, amed, asml, ax, ymax, ymin;

  // Quick return if possible

  if (disnan(scale.value) || disnan(sumsq.value)) return;
  if (sumsq.value == dzero) scale.value = one;
  if (scale.value == zero) {
    scale.value = one;
    sumsq.value = zero;
  }
  if (n <= 0) {
    return;
  }

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
    ax = x[ix].abs();
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
  //
  // Put the existing sum of squares into one of the accumulators
  //
  if (sumsq.value > zero) {
    ax = scale.value * sqrt(sumsq.value);
    if (ax > tbig) {
      if (scale.value > one) {
        scale.value *= sbig;
        abig += scale.value * (scale.value * sumsq.value);
      } else {
        // sumsq > tbig^2 => (sbig * (sbig * sumsq)) is representable;
        abig =
            abig + scale.value * (scale.value * (sbig * (sbig * sumsq.value)));
      }
    } else if (ax < tsml) {
      if (notbig) {
        if (scale.value < one) {
          scale.value *= ssml;
          asml += scale.value * (scale.value * sumsq.value);
        } else {
          // sumsq < tsml^2 => (ssml * (ssml * sumsq)) is representable
          asml += scale.value * (scale.value * (ssml * (ssml * sumsq.value)));
        }
      }
    } else {
      amed += scale.value * (scale.value * sumsq.value);
    }
  }
  //
  // Combine abig and amed or amed and asml if more than one
  // accumulator was used.
  //
  if (abig > zero) {
    //
    //    Combine abig and amed if abig > 0.
    //
    if (amed > zero || amed.isNaN) {
      abig += (amed * sbig) * sbig;
    }
    scale.value = one / sbig;
    sumsq.value = abig;
  } else if (asml > zero) {
    //
    //    Combine amed and asml if asml > 0.
    //
    if (amed > zero || amed.isNaN) {
      amed = sqrt(amed);
      asml = sqrt(asml) / ssml;
      if (asml > amed) {
        ymin = amed;
        ymax = asml;
      } else {
        ymin = asml;
        ymax = amed;
      }
      scale.value = one;
      sumsq.value = pow(ymax, 2) * (one + pow((ymin / ymax), 2));
    } else {
      scale.value = one / ssml;
      sumsq.value = asml;
    }
  } else {
    //
    //    Otherwise all values are mid-range or zero
    //
    scale.value = one;
    sumsq.value = amed;
  }
}
