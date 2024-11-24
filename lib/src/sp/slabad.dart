// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

void slabad(SMALL, LARGE) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  double LARGE, SMALL;
  // ..

// =====================================================================

  // .. Intrinsic Functions ..
  // INTRINSIC LOG10, SQRT
  // ..

  // If it looks like we're on a Cray, take the square root of
  // SMALL and LARGE to avoid overflow and underflow problems.

  // IF( LOG10( LARGE ) > 2000. ) THEN
  // SMALL = sqrt( SMALL )
  // LARGE = sqrt( LARGE )
  // END IF
}
