// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dptt01(
  final int N,
  final Array<double> D_,
  final Array<double> E_,
  final Array<double> DF_,
  final Array<double> EF_,
  final Array<double> WORK_,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final E = E_.having();
  final DF = DF_.having();
  final EF = EF_.having();
  final WORK = WORK_.having();
  const ONE = 1.0, ZERO = 0.0;

  // Quick return if possible

  if (N <= 0) {
    RESID.value = ZERO;
    return;
  }

  final EPS = dlamch('Epsilon');

  // Construct the difference L*D*L' - A.

  WORK[1] = DF[1] - D[1];
  for (var I = 1; I <= N - 1; I++) {
    final DE = DF[I] * EF[I];
    WORK[N + I] = DE - E[I];
    WORK[1 + I] = DE * EF[I] + DF[I + 1] - D[I + 1];
  }

  // Compute the 1-norms of the tridiagonal matrices A and WORK.

  double ANORM;
  if (N == 1) {
    ANORM = D[1];
    RESID.value = WORK[1].abs();
  } else {
    ANORM = max(D[1] + E[1].abs(), D[N] + E[N - 1].abs());
    RESID.value = max(WORK[1].abs() + WORK[N + 1].abs(),
        WORK[N].abs() + WORK[2 * N - 1].abs());
    for (var I = 2; I <= N - 1; I++) {
      ANORM = max(ANORM, D[I] + E[I].abs() + E[I - 1].abs());
      RESID.value = max(RESID.value,
          WORK[I].abs() + WORK[N + I - 1].abs() + WORK[N + I].abs());
    }
  }

  // Compute norm(L*D*L' - A) / (n * norm(A) * EPS)

  if (ANORM <= ZERO) {
    if (RESID.value != ZERO) RESID.value = ONE / EPS;
  } else {
    RESID.value = ((RESID.value / N) / ANORM) / EPS;
  }
}
