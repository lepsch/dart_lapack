// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

REAL scabs1(Z) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  Complex Z;
  // ..

// =====================================================================

  // .. Intrinsic Functions ..
  // INTRINSIC ABS,AIMAG,REAL
  // ..
  SCABS1 = (REAL(Z)).abs() + (AIMAG(Z)).abs();
}
