// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

/* origin: FreeBSD /usr/src/lib/msun/src/e_log10.c */
/*
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunSoft, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */
/*
 * Return the base 10 logarithm of x.  See log.c for most comments.
 *
 * Reduce x to 2^k (1+f) and calculate r = log(1+f) - f + f*f/2
 * as in log.c, then combine and scale in extra precision:
 *    log10(x) = (f - f*f/2 + r)/log(10) + k*log10(2)
 */

import 'dart:typed_data';

const _ivln10hi = 4.34294481878168880939e-01; /* 0x3fdbcb7b, 0x15200000 */
const _ivln10lo = 2.50829467116452752298e-11; /* 0x3dbb9438, 0xca9aadd5 */
const _log10_2hi = 3.01029995663611771306e-01; /* 0x3FD34413, 0x509F6000 */
const _log10_2lo = 3.69423907715893078616e-13; /* 0x3D59FEF3, 0x11F12B36 */
const _Lg1 = 6.666666666666735130e-01; /* 3FE55555 55555593 */
const _Lg2 = 3.999999999940941908e-01; /* 3FD99999 9997FA04 */
const _Lg3 = 2.857142874366239149e-01; /* 3FD24924 94229359 */
const _Lg4 = 2.222219843214978396e-01; /* 3FCC71C5 1D8E78AF */
const _Lg5 = 1.818357216161805012e-01; /* 3FC74664 96CB03DE */
const _Lg6 = 1.531383769920937332e-01; /* 3FC39A09 D078C69F */
const _Lg7 = 1.479819860511658591e-01; /* 3FC2F112 DF3E5244 */
const _0x1p54 = [0x43, 0x50, 00, 00, 00, 00, 00, 00];

class _Union {
  final ByteBuffer _buf;
  _Union(double x) : _buf = Float64List.fromList([x]).buffer;

  double get f => _buf.asFloat64List(0).first;
  set f(double v) => _buf.asFloat64List(0).first = v;

  int get i => _buf.asUint64List(0).first;
  set i(int v) => _buf.asUint64List(0).first = v;
}

double log10(double x) {
  final u = _Union(x);
  double hfsq, f, s, z, R, w, t1, t2, dk, y, hi, lo, val_hi, val_lo;
  int hx; // unsigned
  int k;

  hx = u.i >>> 32;
  k = 0;
  if (hx < 0x00100000 || hx >>> 31 != 0) {
    if (u.i << 1 == 0) {
      return -1 / (x * x); /* log(+-0)=-inf */
    }
    if (hx >>> 31 != 0) {
      return (x - x) / 0.0; /* log(-#) = NaN */
    }
    /* subnormal number, scale x up */
    k -= 54;
    x *= _0x1p54.asDouble();
    u.f = x;
    hx = u.i >>> 32;
  } else if (hx >= 0x7ff00000) {
    return x;
  } else if (hx == 0x3ff00000 && u.i << 32 == 0) {
    return 0;
  }

  /* reduce x into [sqrt(2)/2, sqrt(2)] */
  hx += 0x3ff00000 - 0x3fe6a09e;
  k += (hx >>> 20).toInt() - 0x3ff;
  hx = (hx & 0x000fffff) + 0x3fe6a09e;
  u.i = (hx << 32).toUnsigned(64) | (u.i & 0xffffffff);
  x = u.f;

  f = x - 1.0;
  hfsq = 0.5 * f * f;
  s = f / (2.0 + f);
  z = s * s;
  w = z * z;
  t1 = w * (_Lg2 + w * (_Lg4 + w * _Lg6));
  t2 = z * (_Lg1 + w * (_Lg3 + w * (_Lg5 + w * _Lg7)));
  R = t2 + t1;

  /* See log2.c for details. */
  /* hi+lo = f - hfsq + s*(hfsq+R) ~ log(1+f) */
  hi = f - hfsq;
  u.f = hi;
  u.i &= (-1 << 32).toUnsigned(64);
  hi = u.f;
  lo = f - hi - hfsq + s * (hfsq + R);

  /* val_hi+val_lo ~ log10(1+f) + k*log10(2) */
  val_hi = hi * _ivln10hi;
  dk = k.toDouble();
  y = dk * _log10_2hi;
  val_lo = dk * _log10_2lo + (lo + hi) * _ivln10lo + lo * _ivln10hi;

  /*
   * Extra precision in for adding y is not strictly needed
   * since there is no very large cancellation near x = sqrt(2) or
   * x = 1/sqrt(2), but we do it anyway since it costs little on CPUs
   * with some parallelism and it reduces the error for many args.
   */
  w = y + val_hi;
  val_lo += (y - w) + val_hi;
  val_hi = w;

  return val_lo + val_hi;
}

extension on List<int> {
  double asDouble() =>
      Uint32List.fromList(reversed.toList()).buffer.asFloat64List(0).first;
}
