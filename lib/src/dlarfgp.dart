// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/dnrm2.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlapy2.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/intrinsics/sign.dart';
import 'package:dart_lapack/src/matrix.dart';

void dlarfgp(
  final int N,
  final Box<double> ALPHA,
  final Array<double> X_,
  final int INCX,
  final Box<double> TAU,
) {
  final X = X_.having();
  const TWO = 2.0, ONE = 1.0, ZERO = 0.0;
  int J, KNT = 0;
  double BETA = 0, BIGNUM = 0, EPS, SAVEALPHA, SMLNUM = 0, XNORM;

  if (N <= 0) {
    TAU.value = ZERO;
    return;
  }

  EPS = dlamch('Precision');
  XNORM = dnrm2(N - 1, X, INCX);

  if (XNORM <= EPS * ALPHA.value.abs()) {
    // H  =  [+/-1, 0; I], sign chosen so ALPHA >= 0.
    if (ALPHA.value >= ZERO) {
      // When TAU == ZERO, the vector is special-cased to be
      // all zeros in the application routines.  We do not need
      // to clear it.
      TAU.value = ZERO;
    } else {
      // However, the application routines rely on explicit
      // zero checks when TAU != ZERO, and we must clear X.
      TAU.value = TWO;
      for (J = 1; J <= N - 1; J++) {
        X[1 + (J - 1) * INCX] = 0;
      }
      ALPHA.value = -ALPHA.value;
    }
  } else {
    // general case
    BETA = sign(dlapy2(ALPHA.value, XNORM), ALPHA.value);
    SMLNUM = dlamch('S') / dlamch('E');
    KNT = 0;
    if (BETA.abs() < SMLNUM) {
      // XNORM, BETA may be inaccurate; scale X and recompute them
      BIGNUM = ONE / SMLNUM;
      do {
        KNT++;
        dscal(N - 1, BIGNUM, X, INCX);
        BETA *= BIGNUM;
        ALPHA.value *= BIGNUM;
      } while ((BETA.abs() < SMLNUM) && (KNT < 20));

      // New BETA is at most 1, at least SMLNUM
      XNORM = dnrm2(N - 1, X, INCX);
      BETA = sign(dlapy2(ALPHA.value, XNORM), ALPHA.value);
    }
    SAVEALPHA = ALPHA.value;
    ALPHA.value += BETA;
    if (BETA < ZERO) {
      BETA = -BETA;
      TAU.value = -ALPHA.value / BETA;
    } else {
      ALPHA.value = XNORM * (XNORM / ALPHA.value);
      TAU.value = ALPHA.value / BETA;
      ALPHA.value = -ALPHA.value;
    }

    if (TAU.value.abs() <= SMLNUM) {
      // In the case where the computed TAU ends up being a denormalized number,
      // it loses relative accuracy. This is a BIG problem. Solution: flush TAU
      // to ZERO. This explains the next IF statement.
      if (SAVEALPHA >= ZERO) {
        TAU.value = ZERO;
      } else {
        TAU.value = TWO;
        for (J = 1; J <= N - 1; J++) {
          X[1 + (J - 1) * INCX] = 0;
        }
        BETA = -SAVEALPHA;
      }
    } else {
      // This is the general case.
      dscal(N - 1, ONE / ALPHA.value, X, INCX);
    }

    // If BETA is subnormal, it may lose relative accuracy
    for (J = 1; J <= KNT; J++) {
      BETA *= SMLNUM;
    }
    ALPHA.value = BETA;
  }
}
