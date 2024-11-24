// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void zdrscl(
  final int N,
  final double SA,
  final Array<Complex> SX_,
  final int INCX,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final SX = SX_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool DONE;
  double BIGNUM, CDEN, CDEN1, CNUM, CNUM1, MUL, SMLNUM;

  // Quick return if possible

  if (N <= 0) return;

  // Get machine parameters

  SMLNUM = dlamch('S');
  BIGNUM = ONE / SMLNUM;

  // Initialize the denominator to SA and the numerator to 1.

  CDEN = SA;
  CNUM = ONE;

  do {
    CDEN1 = CDEN * SMLNUM;
    CNUM1 = CNUM / BIGNUM;
    if (CDEN1.abs() > CNUM.abs() && CNUM != ZERO) {
      // Pre-multiply X by SMLNUM if CDEN is large compared to CNUM.

      MUL = SMLNUM;
      DONE = false;
      CDEN = CDEN1;
    } else if (CNUM1.abs() > CDEN.abs()) {
      // Pre-multiply X by BIGNUM if CDEN is small compared to CNUM.

      MUL = BIGNUM;
      DONE = false;
      CNUM = CNUM1;
    } else {
      // Multiply X by CNUM / CDEN and return.

      MUL = CNUM / CDEN;
      DONE = true;
    }

    // Scale the vector X by MUL

    zdscal(N, MUL, SX, INCX);
  } while (!DONE);
}
