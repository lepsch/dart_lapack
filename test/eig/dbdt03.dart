// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dasum.dart';
import 'package:dart_lapack/src/blas/dgemv.dart';
import 'package:dart_lapack/src/blas/idamax.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void dbdt03(
  final String UPLO,
  final int N,
  final int KD,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<double> U_,
  final int LDU,
  final Array<double> S_,
  final Matrix<double> VT_,
  final int LDVT,
  final Array<double> WORK_,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final E = E_.having();
  final U = U_.having(ld: LDU);
  final S = S_.having();
  final VT = VT_.having(ld: LDVT);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  int I, J;
  double BNORM, EPS;

  // Quick return if possible

  RESID.value = ZERO;
  if (N <= 0) return;

  // Compute B - U * S * V' one column at a time.

  BNORM = ZERO;
  if (KD >= 1) {
    // B is bidiagonal.

    if (lsame(UPLO, 'U')) {
      // B is upper bidiagonal.

      for (J = 1; J <= N; J++) {
        for (I = 1; I <= N; I++) {
          WORK[N + I] = S[I] * VT[I][J];
        }
        dgemv(
            'No transpose', N, N, -ONE, U, LDU, WORK(N + 1), 1, ZERO, WORK, 1);
        WORK[J] += D[J];
        if (J > 1) {
          WORK[J - 1] += E[J - 1];
          BNORM = max(BNORM, D[J].abs() + E[J - 1].abs());
        } else {
          BNORM = max(BNORM, D[J].abs());
        }
        RESID.value = max(RESID.value, dasum(N, WORK, 1));
      }
    } else {
      // B is lower bidiagonal.

      for (J = 1; J <= N; J++) {
        for (I = 1; I <= N; I++) {
          WORK[N + I] = S[I] * VT[I][J];
        }
        dgemv(
            'No transpose', N, N, -ONE, U, LDU, WORK(N + 1), 1, ZERO, WORK, 1);
        WORK[J] += D[J];
        if (J < N) {
          WORK[J + 1] += E[J];
          BNORM = max(BNORM, D[J].abs() + E[J].abs());
        } else {
          BNORM = max(BNORM, D[J].abs());
        }
        RESID.value = max(RESID.value, dasum(N, WORK, 1));
      }
    }
  } else {
    // B is diagonal.

    for (J = 1; J <= N; J++) {
      for (I = 1; I <= N; I++) {
        WORK[N + I] = S[I] * VT[I][J];
      }
      dgemv('No transpose', N, N, -ONE, U, LDU, WORK(N + 1), 1, ZERO, WORK, 1);
      WORK[J] += D[J];
      RESID.value = max(RESID.value, dasum(N, WORK, 1));
    }
    J = idamax(N, D, 1);
    BNORM = D[J].abs();
  }

  // Compute norm(B - U * S * V') / ( n * norm(B) * EPS )

  EPS = dlamch('Precision');

  if (BNORM <= ZERO) {
    if (RESID.value != ZERO) RESID.value = ONE / EPS;
  } else {
    if (BNORM >= RESID.value) {
      RESID.value = (RESID.value / BNORM) / (N * EPS);
    } else {
      if (BNORM < ONE) {
        RESID.value = (min(RESID.value, N * BNORM) / BNORM) / (N * EPS);
      } else {
        RESID.value = min(RESID.value / BNORM, N) / (N * EPS);
      }
    }
  }
}
