// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/izamax.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void zget04(
  final int N,
  final int NRHS,
  final Matrix<Complex> X_,
  final int LDX,
  final Matrix<Complex> XACT_,
  final int LDXACT,
  final double RCOND,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having(ld: LDX);
  final XACT = XACT_.having(ld: LDXACT);
  const ZERO = 0.0;

  // Quick exit if N = 0 or NRHS = 0.

  if (N <= 0 || NRHS <= 0) {
    RESID.value = ZERO;
    return;
  }

  // Exit with RESID = 1/EPS if RCOND is invalid.

  final EPS = dlamch('Epsilon');
  if (RCOND < ZERO) {
    RESID.value = 1.0 / EPS;
    return;
  }

  // Compute the maximum of
  //    norm(X - XACT) / ( norm(XACT) * EPS )
  // over all the vectors X and XACT .

  RESID.value = ZERO;
  for (var J = 1; J <= NRHS; J++) {
    final IX = izamax(N, XACT(1, J).asArray(), 1);
    final XNORM = XACT[IX][J].cabs1();
    var DIFFNM = ZERO;
    for (var I = 1; I <= N; I++) {
      DIFFNM = max(DIFFNM, (X[I][J] - XACT[I][J]).cabs1());
    }
    if (XNORM <= ZERO) {
      if (DIFFNM > ZERO) RESID.value = 1.0 / EPS;
    } else {
      RESID.value = max(RESID.value, (DIFFNM / XNORM) * RCOND);
    }
  }
  if (RESID.value * EPS < 1.0) RESID.value /= EPS;
}
