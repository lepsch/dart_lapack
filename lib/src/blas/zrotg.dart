// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/intrinsics/maxexponent.dart';
import 'package:dart_lapack/src/intrinsics/minexponent.dart';
import 'package:dart_lapack/src/intrinsics/radix.dart';

void zrotg(
  final Box<Complex> a,
  final Complex b,
  final Box<double> c,
  final Box<Complex> s,
) {
//  -- Reference BLAS level1 routine --
//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const zero = 0.0, one = 1.0;
  final safmin =
      pow(radix(0.0), max(minexponent(0.0) - 1, 1 - maxexponent(0.0)));
  final safmax =
      pow(radix(0.0), max(1 - minexponent(0.0), maxexponent(0.0) - 1));
  final rtmin = sqrt(safmin);

  Complex r;
  final f = a.value;
  final g = b;
  if (g == Complex.zero) {
    c.value = one;
    s.value = Complex.zero;
    r = f;
  } else if (f == Complex.zero) {
    c.value = zero;
    if (g.real == zero) {
      r = g.imaginary.abs().toComplex();
      s.value = g.conjugate() / r;
    } else if (g.imaginary == zero) {
      r = g.real.abs().toComplex();
      s.value = g.conjugate() / r;
    } else {
      final g1 = max(g.real.abs(), g.imaginary.abs());
      final rtmax = sqrt(safmax / 2);
      if (g1 > rtmin && g1 < rtmax) {
        // Use unscaled algorithm
        //
        //    The following two lines can be replaced by `d = abs( g )`.
        //    This algorithm do not use the intrinsic complex abs.

        final g2 = g.cabsSq();
        final d = sqrt(g2);
        s.value = g.conjugate() / d.toComplex();
        r = d.toComplex();
      } else {
        // Use scaled algorithm

        final u = min(safmax, max(safmin, g1));
        final gs = g / u.toComplex();
        // The following two lines can be replaced by `d = abs( gs )`.
        // This algorithm do not use the intrinsic complex abs.
        final g2 = gs.cabsSq();
        final d = sqrt(g2);
        s.value = gs.conjugate() / d.toComplex();
        r = (d * u).toComplex();
      }
    }
  } else {
    final f1 = max(f.real.abs(), f.imaginary.abs());
    final g1 = max(g.real.abs(), g.imaginary.abs());
    var rtmax = sqrt(safmax / 4);
    if (f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax) {
      // Use unscaled algorithm

      final f2 = f.cabsSq();
      final g2 = g.cabsSq();
      final h2 = f2 + g2;
      // safmin <= f2 <= h2 <= safmax
      if (f2 >= h2 * safmin) {
        // safmin <= f2/h2 <= 1, and h2/f2 is finite
        c.value = sqrt(f2 / h2);
        r = f / c.value.toComplex();
        rtmax *= 2;
        if (f2 > rtmin && h2 < rtmax) {
          // safmin <= sqrt( f2*h2 ) <= safmax
          s.value = g.conjugate() * (f / sqrt(f2 * h2).toComplex());
        } else {
          s.value = g.conjugate() * (r / h2.toComplex());
        }
      } else {
        // f2/h2 <= safmin may be subnormal, and h2/f2 may overflow.
        // Moreover,
        //  safmin <= f2*f2 * safmax < f2 * h2 < h2*h2 * safmin <= safmax,
        //  sqrt(safmin) <= sqrt(f2 * h2) <= sqrt(safmax).
        // Also,
        //  g2 >> f2, which means that h2 = g2.
        final d = sqrt(f2 * h2);
        c.value = f2 / d;
        if (c.value >= safmin) {
          r = f / c.value.toComplex();
        } else {
          // f2 / sqrt(f2 * h2) < safmin, {
          //  sqrt(safmin) <= f2 * sqrt(safmax) <= h2 / sqrt(f2 * h2) <= h2 * (safmin / f2) <= h2 <= safmax
          r = f * (h2 / d).toComplex();
        }
        s.value = g.conjugate() * (f / d.toComplex());
      }
    } else {
      // Use scaled algorithm

      final u = min(safmax, max(safmin, max(f1, g1)));
      final gs = g / u.toComplex();
      final g2 = gs.cabsSq();
      final double w, f2, h2;
      final Complex fs;
      if (f1 / u < rtmin) {
        // f is not well-scaled when scaled by g1.
        // Use a different scaling for f.

        final v = min(safmax, max(safmin, f1));
        w = v / u;
        fs = f / v.toComplex();
        f2 = fs.cabsSq();
        h2 = f2 * pow(w, 2) + g2;
      } else {
        // Otherwise use the same scaling for f and g.

        w = one;
        fs = f / u.toComplex();
        f2 = fs.cabsSq();
        h2 = f2 + g2;
      }
      // safmin <= f2 <= h2 <= safmax
      if (f2 >= h2 * safmin) {
        // safmin <= f2/h2 <= 1, and h2/f2 is finite
        c.value = sqrt(f2 / h2);
        r = fs / c.value.toComplex();
        rtmax *= 2;
        if (f2 > rtmin && h2 < rtmax) {
          // safmin <= sqrt( f2*h2 ) <= safmax
          s.value = gs.conjugate() * (fs / sqrt(f2 * h2).toComplex());
        } else {
          s.value = gs.conjugate() * (r / h2.toComplex());
        }
      } else {
        // f2/h2 <= safmin may be subnormal, and h2/f2 may overflow.
        // Moreover,
        //  safmin <= f2*f2 * safmax < f2 * h2 < h2*h2 * safmin <= safmax,
        //  sqrt(safmin) <= sqrt(f2 * h2) <= sqrt(safmax).
        // Also,
        //  g2 >> f2, which means that h2 = g2.
        final d = sqrt(f2 * h2);
        c.value = f2 / d;
        if (c.value >= safmin) {
          r = fs / c.value.toComplex();
        } else {
          // f2 / sqrt(f2 * h2) < safmin, {
          //  sqrt(safmin) <= f2 * sqrt(safmax) <= h2 / sqrt(f2 * h2) <= h2 * (safmin / f2) <= h2 <= safmax
          r = fs * (h2 / d).toComplex();
        }
        s.value = gs.conjugate() * (fs / d.toComplex());
      }
      // Rescale c and r
      c.value *= w;
      r *= u.toComplex();
    }
  }
  a.value = r;
}
