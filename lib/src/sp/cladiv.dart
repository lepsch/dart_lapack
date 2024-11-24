// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

Complex cladiv(X, Y) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  Complex X, Y;
  // ..

// =====================================================================

  // .. Local Scalars ..
  REAL ZI, ZR;
  // ..
  // .. External Subroutines ..
  // EXTERNAL SLADIV
  // ..
  // .. Intrinsic Functions ..
  // INTRINSIC AIMAG, CMPLX, REAL
  // ..

  sladiv(REAL(X), AIMAG(X), REAL(Y), AIMAG(Y), ZR, ZI);
  CLADIV = CMPLX(ZR, ZI);
}
