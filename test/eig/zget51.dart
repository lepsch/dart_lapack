// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';

void zget51(
  final int ITYPE,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> U_,
  final int LDU,
  final Matrix<Complex> V_,
  final int LDV,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Box<double> RESULT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final U = U_.having(ld: LDU);
  final V = V_.having(ld: LDV);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
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

    ANORM = max(zlange('1', N, N, A, LDA, RWORK), UNFL);

    if (ITYPE == 1) {
      // ITYPE=1: Compute W = A - U B V**H

      zlacpy(' ', N, N, A, LDA, WORK.asMatrix(), N);
      zgemm('N', 'N', N, N, N, Complex.one, U, LDU, B, LDB, Complex.zero,
          WORK(pow(N, 2).toInt() + 1).asMatrix(), N);

      zgemm(
          'N',
          'C',
          N,
          N,
          N,
          -Complex.one,
          WORK(pow(N, 2).toInt() + 1).asMatrix(),
          N,
          V,
          LDV,
          Complex.one,
          WORK.asMatrix(),
          N);
    } else {
      // ITYPE=2: Compute W = A - B

      zlacpy(' ', N, N, B, LDB, WORK.asMatrix(), N);

      for (JCOL = 1; JCOL <= N; JCOL++) {
        for (JROW = 1; JROW <= N; JROW++) {
          WORK[JROW + N * (JCOL - 1)] =
              WORK[JROW + N * (JCOL - 1)] - A[JROW][JCOL];
        }
      }
    }

    // Compute norm(W)/ ( ulp*norm(A) )

    WNORM = zlange('1', N, N, WORK.asMatrix(), N, RWORK);

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

    // ITYPE=3: Compute  U U**H - I

    zgemm('N', 'C', N, N, N, Complex.one, U, LDU, U, LDU, Complex.zero,
        WORK.asMatrix(), N);

    for (JDIAG = 1; JDIAG <= N; JDIAG++) {
      WORK[(N + 1) * (JDIAG - 1) + 1] =
          WORK[(N + 1) * (JDIAG - 1) + 1] - Complex.one;
    }

    RESULT.value =
        min(zlange('1', N, N, WORK.asMatrix(), N, RWORK), N) / (N * ULP);
  }
}
