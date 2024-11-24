// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

double slarmm(ANORM, BNORM, CNORM) {
  // IMPLICIT NONE
  // .. Scalar Arguments ..
  double ANORM, BNORM, CNORM;
  // .. Parameters ..
  double ONE, HALF, FOUR;
  const ONE = 1.0, HALF = 0.5, FOUR = 4.0;
  // ..
  // .. Local Scalars ..
  double BIGNUM, SMLNUM;
  // ..
  // .. External Functions ..
  //- REAL               SLAMCH;
  // EXTERNAL SLAMCH
  // ..

  // Determine machine dependent parameters to control overflow.

  SMLNUM = SLAMCH('Safe minimum') / SLAMCH('Precision');
  BIGNUM = (ONE / SMLNUM) / FOUR;

  // Compute a scale factor.

  SLARMM = ONE;
  if (BNORM <= ONE) {
    if (ANORM * BNORM > BIGNUM - CNORM) {
      SLARMM = HALF;
    }
  } else {
    if (ANORM > (BIGNUM - CNORM) / BNORM) {
      SLARMM = HALF / BNORM;
    }
  }
}
