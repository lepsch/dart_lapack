// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/blas/zhemm.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlanhe.dart';

import 'zunt01.dart';

void zhet22(
  final int ITYPE,
  final String UPLO,
  final int N,
  final int M,
  final int KBAND,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<Complex> U_,
  final int LDU,
  final Matrix<Complex> V_,
  final int LDV,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Array<double> RESULT_,
) {
  final A = A_.having(ld: LDA);
  final U = U_.having(ld: LDU);
  // final V = V_.having(ld: LDV);
  // final TAU = TAU_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final D = D_.having();
  final E = E_.having();
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

  ANORM = max(zlanhe('1', UPLO, N, A, LDA, RWORK), UNFL);

  // Compute error matrix:

  // ITYPE=1: error = U**H A U - S

  zhemm('L', UPLO, N, M, Complex.one, A, LDA, U, LDU, Complex.zero,
      WORK.asMatrix(), N);
  NN = N * N;
  NNP1 = NN + 1;
  zgemm('C', 'N', M, M, N, Complex.one, U, LDU, WORK.asMatrix(), N,
      Complex.zero, WORK(NNP1).asMatrix(), N);
  for (J = 1; J <= M; J++) {
    JJ = NN + (J - 1) * N + J;
    WORK[JJ] -= D[J].toComplex();
  }
  if (KBAND == 1 && N > 1) {
    for (J = 2; J <= M; J++) {
      JJ1 = NN + (J - 1) * N + J - 1;
      JJ2 = NN + (J - 2) * N + J;
      WORK[JJ1] -= E[J - 1].toComplex();
      WORK[JJ2] -= E[J - 1].toComplex();
    }
  }
  WNORM = zlanhe('1', UPLO, M, WORK(NNP1).asMatrix(), N, RWORK);

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

  // Compute  U**H U - I

  if (ITYPE == 1) {
    zunt01('Columns', N, M, U, LDU, WORK, 2 * N * N, RWORK, RESULT(2));
  }
}
