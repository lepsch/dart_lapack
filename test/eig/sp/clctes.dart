// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

bool clctes(Z, D) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  Complex D, Z;
  // ..

// =====================================================================

  // .. Parameters ..

  REAL ZERO, ONE;
  const ZERO = 0.0, ONE = 1.0;
  Complex CZERO;
  const CZERO = (0.0, 0.0);
  // ..
  // .. Local Scalars ..
  REAL ZMAX;
  // ..
  // .. Intrinsic Functions ..
  // INTRINSIC ABS, AIMAG, MAX, REAL, SIGN
  // ..

  if (D == CZERO) {
    CLCTES = (REAL(Z) < ZERO);
  } else {
    if (REAL(Z) == ZERO || REAL(D) == ZERO) {
      CLCTES = (sign(ONE, AIMAG(Z)) != sign(ONE, AIMAG(D)));
    } else if (AIMAG(Z) == ZERO || AIMAG(D) == ZERO) {
      CLCTES = (sign(ONE, REAL(Z)) != sign(ONE, double(D)));
    } else {
      ZMAX = max((double(Z)).abs(), (AIMAG(Z)).abs());
      CLCTES =
          ((double(Z) / ZMAX) * REAL(D) + (AIMAG(Z) / ZMAX) * AIMAG(D) < ZERO);
    }
  }
}
