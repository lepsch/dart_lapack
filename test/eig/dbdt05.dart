// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dasum.dart';
import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void dbdt05(
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> S_,
  final int NS,
  final Matrix<double> U_,
  final int LDU,
  final Matrix<double> VT_,
  final int LDVT,
  final Array<double> WORK_,
  final Box<double> RESID,
) {
  final A = A_.having(ld: LDA);
  final S = S_.having();
  final U = U_.having(ld: LDU);
  final VT = VT_.having(ld: LDVT);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  int I, J;
  double ANORM, EPS;

  // Quick return if possible.

  RESID.value = ZERO;
  if (min(M, N) <= 0 || NS <= 0) return;

  EPS = dlamch('Precision');
  ANORM = dlange('M', M, N, A, LDA, WORK);

  // Compute U' * A * V.

  dgemm('N', 'T', M, NS, N, ONE, A, LDA, VT, LDVT, ZERO,
      WORK(1 + NS * NS).asMatrix(M), M);
  dgemm('T', 'N', NS, NS, M, -ONE, U, LDU, WORK(1 + NS * NS).asMatrix(M), M,
      ZERO, WORK.asMatrix(NS), NS);

  // norm(S - U' * B * V)

  J = 0;
  for (I = 1; I <= NS; I++) {
    WORK[J + I] += S[I];
    RESID.value = max(RESID.value, dasum(NS, WORK(J + 1), 1));
    J += NS;
  }

  if (ANORM <= ZERO) {
    if (RESID.value != ZERO) RESID.value = ONE / EPS;
  } else {
    if (ANORM >= RESID.value) {
      RESID.value = (RESID.value / ANORM) / (N * EPS);
    } else {
      if (ANORM < ONE) {
        RESID.value = (min(RESID.value, N * ANORM) / ANORM) / (N * EPS);
      } else {
        RESID.value = min(RESID.value / ANORM, N) / (N * EPS);
      }
    }
  }
}
