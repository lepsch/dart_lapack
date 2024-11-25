// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlantr.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void dtrt06(
  final double RCOND,
  final double RCONDC,
  final String UPLO,
  final String DIAG,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> WORK_,
  final Box<double> RAT,
) {
  final A = A_.having(ld: LDA);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  final EPS = dlamch('Epsilon');
  final RMAX = max(RCOND, RCONDC);
  final RMIN = min(RCOND, RCONDC);

  // Do the easy cases first.

  if (RMIN < ZERO) {
    // Invalid value for RCOND or RCONDC, return 1/EPS.

    RAT.value = ONE / EPS;
  } else if (RMIN > ZERO) {
    // Both estimates are positive, return RMAX/RMIN - 1.

    RAT.value = RMAX / RMIN - ONE;
  } else if (RMAX == ZERO) {
    // Both estimates zero.

    RAT.value = ZERO;
  } else {
    // One estimate is zero, the other is non-zero.  If the matrix is
    // ill-conditioned, return the nonzero estimate multiplied by
    // 1/EPS; if the matrix is badly scaled, return the nonzero
    // estimate multiplied by BIGNUM/TMAX, where TMAX is the maximum
    // element in absolute value in A.

    final SMLNUM = dlamch('Safe minimum');
    final BIGNUM = ONE / SMLNUM;
    final ANORM = dlantr('M', UPLO, DIAG, N, N, A, LDA, WORK);

    RAT.value = RMAX * (min(BIGNUM / max(ONE, ANORM), ONE / EPS));
  }
}
