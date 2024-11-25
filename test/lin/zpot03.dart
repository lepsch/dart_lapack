// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/zhemm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlanhe.dart';

void zpot03(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> AINV_,
  final int LDAINV,
  final Matrix<Complex> WORK_,
  final int LDWORK,
  final Array<double> RWORK_,
  final Box<double> RCOND,
  final Box<double> RESID,
) {
  final A = A_.having(ld: LDA);
  final AINV = AINV_.having(ld: LDAINV);
  final WORK = WORK_.having(ld: LDWORK);
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if N = 0.

  if (N <= 0) {
    RCOND.value = ONE;
    RESID.value = ZERO;
    return;
  }

  // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

  final EPS = dlamch('Epsilon');
  final ANORM = zlanhe('1', UPLO, N, A, LDA, RWORK);
  final AINVNM = zlanhe('1', UPLO, N, AINV, LDAINV, RWORK);
  if (ANORM <= ZERO || AINVNM <= ZERO) {
    RCOND.value = ZERO;
    RESID.value = ONE / EPS;
    return;
  }
  RCOND.value = (ONE / ANORM) / AINVNM;

  // Expand AINV into a full matrix and call ZHEMM to multiply
  // AINV on the left by A.

  if (lsame(UPLO, 'U')) {
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= J - 1; I++) {
        AINV[J][I] = AINV[I][J].conjugate();
      }
    }
  } else {
    for (var J = 1; J <= N; J++) {
      for (var I = J + 1; I <= N; I++) {
        AINV[J][I] = AINV[I][J].conjugate();
      }
    }
  }
  zhemm('Left', UPLO, N, N, -Complex.one, A, LDA, AINV, LDAINV, Complex.zero,
      WORK, LDWORK);

  // Add the identity matrix to WORK .

  for (var I = 1; I <= N; I++) {
    WORK[I][I] += Complex.one;
  }

  // Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)

  RESID.value = zlange('1', N, N, WORK, LDWORK, RWORK);

  RESID.value = ((RESID.value * RCOND.value) / EPS) / N;
}
