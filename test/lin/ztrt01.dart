// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/ztrmv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlantr.dart';

void ztrt01(
  final String UPLO,
  final String DIAG,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> AINV_,
  final int LDAINV,
  final Box<double> RCOND,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
  final A = A_.having(ld: LDA);
  final AINV = AINV_.having(ld: LDAINV);
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if N = 0

  if (N <= 0) {
    RCOND.value = ONE;
    RESID.value = ZERO;
    return;
  }

  // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

  final EPS = dlamch('Epsilon');
  final ANORM = zlantr('1', UPLO, DIAG, N, N, A, LDA, RWORK);
  final AINVNM = zlantr('1', UPLO, DIAG, N, N, AINV, LDAINV, RWORK);
  if (ANORM <= ZERO || AINVNM <= ZERO) {
    RCOND.value = ZERO;
    RESID.value = ONE / EPS;
    return;
  }
  RCOND.value = (ONE / ANORM) / AINVNM;

  // Set the diagonal of AINV to 1 if AINV has unit diagonal.

  if (lsame(DIAG, 'U')) {
    for (var J = 1; J <= N; J++) {
      AINV[J][J] = Complex.one;
    }
  }

  // Compute A * AINV, overwriting AINV.

  if (lsame(UPLO, 'U')) {
    for (var J = 1; J <= N; J++) {
      ztrmv('Upper', 'No transpose', DIAG, J, A, LDA, AINV(1, J).asArray(), 1);
    }
  } else {
    for (var J = 1; J <= N; J++) {
      ztrmv('Lower', 'No transpose', DIAG, N - J + 1, A(J, J), LDA,
          AINV(J, J).asArray(), 1);
    }
  }

  // Subtract 1 from each diagonal element to form A*AINV - I.

  for (var J = 1; J <= N; J++) {
    AINV[J][J] -= Complex.one;
  }

  // Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)

  RESID.value = zlantr('1', UPLO, 'Non-unit', N, N, AINV, LDAINV, RWORK);

  RESID.value = ((RESID.value * RCOND.value) / N) / EPS;
}
