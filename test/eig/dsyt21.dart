// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/dsyr.dart';
import 'package:dart_lapack/src/blas/dsyr2.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dlansy.dart';
import 'package:dart_lapack/src/dlarfy.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dorm2l.dart';
import 'package:dart_lapack/src/dorm2r.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void dsyt21(
  final int ITYPE,
  final String UPLO,
  final int N,
  final int KBAND,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<double> U_,
  final int LDU,
  final Matrix<double> V_,
  final int LDV,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final Array<double> RESULT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final D = D_.having();
  final E = E_.having();
  final U = U_.having(ld: LDU);
  final V = V_.having(ld: LDV);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0, TEN = 10.0;
  bool LOWER;
  String CUPLO;
  int J, JCOL, JR, JROW;
  double ANORM, ULP, UNFL, VSAVE, WNORM = 0;
  final IINFO = Box(0);

  RESULT[1] = ZERO;
  if (ITYPE == 1) RESULT[2] = ZERO;
  if (N <= 0) return;

  if (lsame(UPLO, 'U')) {
    LOWER = false;
    CUPLO = 'U';
  } else {
    LOWER = true;
    CUPLO = 'L';
  }

  UNFL = dlamch('Safe minimum');
  ULP = dlamch('Epsilon') * dlamch('Base');

  // Some Error Checks

  if (ITYPE < 1 || ITYPE > 3) {
    RESULT[1] = TEN / ULP;
    return;
  }

  // Do Test 1

  // Norm of A:

  if (ITYPE == 3) {
    ANORM = ONE;
  } else {
    ANORM = max(dlansy('1', CUPLO, N, A, LDA, WORK), UNFL);
  }

  // Compute error matrix:

  if (ITYPE == 1) {
    // ITYPE=1: error = A - U S U**T

    dlaset('Full', N, N, ZERO, ZERO, WORK.asMatrix(N), N);
    dlacpy(CUPLO, N, N, A, LDA, WORK.asMatrix(N), N);

    for (J = 1; J <= N; J++) {
      dsyr(CUPLO, N, -D[J], U(1, J).asArray(), 1, WORK.asMatrix(N), N);
    }

    if (N > 1 && KBAND == 1) {
      for (J = 1; J <= N - 1; J++) {
        dsyr2(CUPLO, N, -E[J], U(1, J).asArray(), 1, U(1, J + 1).asArray(), 1,
            WORK.asMatrix(N), N);
      }
    }
    WNORM =
        dlansy('1', CUPLO, N, WORK.asMatrix(N), N, WORK(pow(N, 2).toInt() + 1));
  } else if (ITYPE == 2) {
    // ITYPE=2: error = V S V**T - A

    dlaset('Full', N, N, ZERO, ZERO, WORK.asMatrix(N), N);

    if (LOWER) {
      WORK[pow(N, 2).toInt()] = D[N];
      for (J = N - 1; J >= 1; J--) {
        if (KBAND == 1) {
          WORK[(N + 1) * (J - 1) + 2] = (ONE - TAU[J]) * E[J];
          for (JR = J + 2; JR <= N; JR++) {
            WORK[(J - 1) * N + JR] = -TAU[J] * E[J] * V[JR][J];
          }
        }

        VSAVE = V[J + 1][J];
        V[J + 1][J] = ONE;
        dlarfy('L', N - J, V(J + 1, J).asArray(), 1, TAU[J],
            WORK((N + 1) * J + 1).asMatrix(N), N, WORK(pow(N, 2).toInt() + 1));
        V[J + 1][J] = VSAVE;
        WORK[(N + 1) * (J - 1) + 1] = D[J];
      }
    } else {
      WORK[1] = D[1];
      for (J = 1; J <= N - 1; J++) {
        if (KBAND == 1) {
          WORK[(N + 1) * J] = (ONE - TAU[J]) * E[J];
          for (JR = 1; JR <= J - 1; JR++) {
            WORK[J * N + JR] = -TAU[J] * E[J] * V[JR][J + 1];
          }
        }

        VSAVE = V[J][J + 1];
        V[J][J + 1] = ONE;
        dlarfy('U', J, V(1, J + 1).asArray(), 1, TAU[J], WORK.asMatrix(N), N,
            WORK(pow(N, 2).toInt() + 1));
        V[J][J + 1] = VSAVE;
        WORK[(N + 1) * J + 1] = D[J + 1];
      }
    }

    for (JCOL = 1; JCOL <= N; JCOL++) {
      if (LOWER) {
        for (JROW = JCOL; JROW <= N; JROW++) {
          WORK[JROW + N * (JCOL - 1)] =
              WORK[JROW + N * (JCOL - 1)] - A[JROW][JCOL];
        }
      } else {
        for (JROW = 1; JROW <= JCOL; JROW++) {
          WORK[JROW + N * (JCOL - 1)] =
              WORK[JROW + N * (JCOL - 1)] - A[JROW][JCOL];
        }
      }
    }
    WNORM =
        dlansy('1', CUPLO, N, WORK.asMatrix(N), N, WORK(pow(N, 2).toInt() + 1));
  } else if (ITYPE == 3) {
    // ITYPE=3: error = U V**T - I

    if (N < 2) return;
    dlacpy(' ', N, N, U, LDU, WORK.asMatrix(N), N);
    if (LOWER) {
      dorm2r('R', 'T', N, N - 1, N - 1, V(2, 1), LDV, TAU,
          WORK(N + 1).asMatrix(N), N, WORK(pow(N, 2).toInt() + 1), IINFO);
    } else {
      dorm2l('R', 'T', N, N - 1, N - 1, V(1, 2), LDV, TAU, WORK.asMatrix(N), N,
          WORK(pow(N, 2).toInt() + 1), IINFO);
    }
    if (IINFO.value != 0) {
      RESULT[1] = TEN / ULP;
      return;
    }

    for (J = 1; J <= N; J++) {
      WORK[(N + 1) * (J - 1) + 1] -= ONE;
    }

    WNORM = dlange('1', N, N, WORK.asMatrix(N), N, WORK(pow(N, 2).toInt() + 1));
  }

  if (ANORM > WNORM) {
    RESULT[1] = (WNORM / ANORM) / (N * ULP);
  } else {
    if (ANORM < ONE) {
      RESULT[1] = (min(WNORM, N * ANORM) / ANORM) / (N * ULP);
    } else {
      RESULT[1] = min(WNORM / ANORM, N) / (N * ULP);
    }
  }

  // Do Test 2

  // Compute  U U**T - I

  if (ITYPE == 1) {
    dgemm('N', 'C', N, N, N, ONE, U, LDU, U, LDU, ZERO, WORK.asMatrix(N), N);

    for (J = 1; J <= N; J++) {
      WORK[(N + 1) * (J - 1) + 1] -= ONE;
    }

    RESULT[2] = min(
            dlange('1', N, N, WORK.asMatrix(N), N, WORK(pow(N, 2).toInt() + 1)),
            N) /
        (N * ULP);
  }
}
