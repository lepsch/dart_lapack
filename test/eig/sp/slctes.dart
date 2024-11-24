// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

bool slctes(ZR, ZI, D) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  double D, ZI, ZR;
  // ..

// =====================================================================

  // .. Parameters ..
  double ZERO, ONE;
  const ZERO = 0.0, ONE = 1.0;
  // ..
  // .. Intrinsic Functions ..
  // INTRINSIC SIGN
  // ..

  if (D == ZERO) {
    SLCTES = (ZR < ZERO);
  } else {
    SLCTES = (sign(ONE, ZR) != sign(ONE, D));
  }
}
