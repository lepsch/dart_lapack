// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/dsymm.dart';
import 'package:dart_lapack/src/dlansy.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

import 'dort01.dart';

void dsyt22(
  final int ITYPE,
  final String UPLO,
  final int N,
  final int M,
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
  final Array<double> RESULT_,
) {
  final A = A_.having(ld: LDA);
  final D = D_.having();
  final E = E_.having();
  final U = U_.having(ld: LDU);
  // final V = V_.having(ld: LDV);
  // final TAU = TAU_.having();
  final WORK = WORK_.having();
  final RESULT = RESULT_.having();
  const ZERO = 0.0, ONE = 1.0;
  int J, JJ, JJ1, JJ2, NN, NNP1;
  double ANORM, ULP, UNFL, WNORM;

  RESULT[1] = ZERO;
  RESULT[2] = ZERO;
  if (N <= 0 || M <= 0) return;

  UNFL = dlamch('Safe minimum');
  ULP = dlamch('Precision');

  // Do Test 1

  // Norm of A:

  ANORM = max(dlansy('1', UPLO, N, A, LDA, WORK), UNFL);

  // Compute error matrix:

  // ITYPE=1: error = U**T A U - S

  dsymm('L', UPLO, N, M, ONE, A, LDA, U, LDU, ZERO, WORK.asMatrix(N), N);
  NN = N * N;
  NNP1 = NN + 1;
  dgemm('T', 'N', M, M, N, ONE, U, LDU, WORK.asMatrix(N), N, ZERO,
      WORK(NNP1).asMatrix(N), N);
  for (J = 1; J <= M; J++) {
    JJ = NN + (J - 1) * N + J;
    WORK[JJ] -= D[J];
  }
  if (KBAND == 1 && N > 1) {
    for (J = 2; J <= M; J++) {
      JJ1 = NN + (J - 1) * N + J - 1;
      JJ2 = NN + (J - 2) * N + J;
      WORK[JJ1] -= E[J - 1];
      WORK[JJ2] -= E[J - 1];
    }
  }
  WNORM = dlansy('1', UPLO, M, WORK(NNP1).asMatrix(N), N, WORK(1));

  if (ANORM > WNORM) {
    RESULT[1] = (WNORM / ANORM) / (M * ULP);
  } else {
    if (ANORM < ONE) {
      RESULT[1] = (min(WNORM, M * ANORM) / ANORM) / (M * ULP);
    } else {
      RESULT[1] = min(WNORM / ANORM, M) / (M * ULP);
    }
  }

  // Do Test 2

  // Compute  U**T U - I

  if (ITYPE == 1) {
    dort01('Columns', N, M, U, LDU, WORK, 2 * N * N, RESULT.box(2));
  }
}
