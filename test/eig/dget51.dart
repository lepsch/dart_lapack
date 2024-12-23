// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void dget51(
  final int ITYPE,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> U_,
  final int LDU,
  final Matrix<double> V_,
  final int LDV,
  final Array<double> WORK_,
  final Box<double> RESULT,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final U = U_.having(ld: LDU);
  final V = V_.having(ld: LDV);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0, TEN = 10.0;
  int JCOL, JDIAG, JROW;
  double ANORM, ULP, UNFL, WNORM;

  RESULT.value = ZERO;
  if (N <= 0) return;

  // Constants

  UNFL = dlamch('Safe minimum');
  ULP = dlamch('Epsilon') * dlamch('Base');

  // Some Error Checks

  if (ITYPE < 1 || ITYPE > 3) {
    RESULT.value = TEN / ULP;
    return;
  }

  if (ITYPE <= 2) {
    // Tests scaled by the norm(A)

    ANORM = max(dlange('1', N, N, A, LDA, WORK), UNFL);

    if (ITYPE == 1) {
      // ITYPE=1: Compute W = A - UBV'

      dlacpy(' ', N, N, A, LDA, WORK.asMatrix(N), N);
      dgemm('N', 'N', N, N, N, ONE, U, LDU, B, LDB, ZERO,
          WORK(pow(N, 2).toInt() + 1).asMatrix(N), N);

      dgemm('N', 'C', N, N, N, -ONE, WORK(pow(N, 2).toInt() + 1).asMatrix(N), N,
          V, LDV, ONE, WORK.asMatrix(N), N);
    } else {
      // ITYPE=2: Compute W = A - B

      dlacpy(' ', N, N, B, LDB, WORK.asMatrix(N), N);

      for (JCOL = 1; JCOL <= N; JCOL++) {
        for (JROW = 1; JROW <= N; JROW++) {
          WORK[JROW + N * (JCOL - 1)] =
              WORK[JROW + N * (JCOL - 1)] - A[JROW][JCOL];
        }
      }
    }

    // Compute norm(W)/ ( ulp*norm(A) )

    WNORM = dlange('1', N, N, WORK.asMatrix(N), N, WORK(pow(N, 2).toInt() + 1));

    if (ANORM > WNORM) {
      RESULT.value = (WNORM / ANORM) / (N * ULP);
    } else {
      if (ANORM < ONE) {
        RESULT.value = (min(WNORM, N * ANORM) / ANORM) / (N * ULP);
      } else {
        RESULT.value = min(WNORM / ANORM, N) / (N * ULP);
      }
    }
  } else {
    // Tests not scaled by norm(A)

    // ITYPE=3: Compute  UU' - I

    dgemm('N', 'C', N, N, N, ONE, U, LDU, U, LDU, ZERO, WORK.asMatrix(N), N);

    for (JDIAG = 1; JDIAG <= N; JDIAG++) {
      WORK[(N + 1) * (JDIAG - 1) + 1] -= ONE;
    }

    RESULT.value = min(
          dlange('1', N, N, WORK.asMatrix(N), N, WORK(pow(N, 2).toInt() + 1)),
          N,
        ) /
        (N * ULP);
  }
}
