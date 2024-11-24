// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dzasum.dart';
import 'package:dart_lapack/src/blas/zcopy.dart';
import 'package:dart_lapack/src/blas/zgemv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlange.dart';

void zbdt02(
  final int M,
  final int N,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> C_,
  final int LDC,
  final Matrix<Complex> U_,
  final int LDU,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final B = B_.having(ld: LDB);
  final C = C_.having(ld: LDC);
  final U = U_.having(ld: LDU);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  int J;
  double BNORM, EPS, REALMN;

  // Quick return if possible

  RESID.value = ZERO;
  if (M <= 0 || N <= 0) return;
  REALMN = max(M, N).toDouble();
  EPS = dlamch('Precision');

  // Compute norm(B - U * C)

  for (J = 1; J <= N; J++) {
    zcopy(M, B(1, J).asArray(), 1, WORK, 1);
    zgemv('No transpose', M, M, -Complex.one, U, LDU, C(1, J).asArray(), 1,
        Complex.one, WORK, 1);
    RESID.value = max(RESID.value, dzasum(M, WORK, 1));
  }

  // Compute norm of B.

  BNORM = zlange('1', M, N, B, LDB, RWORK);

  if (BNORM <= ZERO) {
    if (RESID.value != ZERO) RESID.value = ONE / EPS;
  } else {
    if (BNORM >= RESID.value) {
      RESID.value = (RESID.value / BNORM) / (REALMN * EPS);
    } else {
      if (BNORM < ONE) {
        RESID.value =
            (min(RESID.value, REALMN * BNORM) / BNORM) / (REALMN * EPS);
      } else {
        RESID.value = min(RESID.value / BNORM, REALMN) / (REALMN * EPS);
      }
    }
  }
}
