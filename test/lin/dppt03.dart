// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dspmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlansp.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';

void dppt03(
  final String UPLO,
  final int N,
  final Array<double> A_,
  final Array<double> AINV_,
  final Matrix<double> WORK_,
  final int LDWORK,
  final Array<double> RWORK_,
  final Box<double> RCOND,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having();
  final AINV = AINV_.having();
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
  final ANORM = dlansp('1', UPLO, N, A, RWORK);
  final AINVNM = dlansp('1', UPLO, N, AINV, RWORK);
  if (ANORM <= ZERO || AINVNM == ZERO) {
    RCOND.value = ZERO;
    RESID.value = ONE / EPS;
    return;
  }
  RCOND.value = (ONE / ANORM) / AINVNM;

  // UPLO = 'U':
  // Copy the leading N-1 x N-1 submatrix of AINV to WORK(1:N,2:N) and
  // expand it to a full matrix, then multiply by A one column at a
  // time, moving the result one column to the left.

  if (lsame(UPLO, 'U')) {
    // Copy AINV

    var JJ = 1;
    for (var J = 1; J <= N - 1; J++) {
      dcopy(J, AINV(JJ), 1, WORK(1, J + 1).asArray(), 1);
      dcopy(J - 1, AINV(JJ), 1, WORK(J, 2).asArray(), LDWORK);
      JJ += J;
    }
    JJ = ((N - 1) * N) ~/ 2 + 1;
    dcopy(N - 1, AINV(JJ), 1, WORK(N, 2).asArray(), LDWORK);

    // Multiply by A

    for (var J = 1; J <= N - 1; J++) {
      dspmv('Upper', N, -ONE, A, WORK(1, J + 1).asArray(), 1, ZERO,
          WORK(1, J).asArray(), 1);
    }
    dspmv('Upper', N, -ONE, A, AINV(JJ), 1, ZERO, WORK(1, N).asArray(), 1);

    // UPLO = 'L':
    // Copy the trailing N-1 x N-1 submatrix of AINV to WORK(1:N,1:N-1)
    // and multiply by A, moving each column to the right.
  } else {
    // Copy AINV

    dcopy(N - 1, AINV(2), 1, WORK(1, 1).asArray(), LDWORK);
    var JJ = N + 1;
    for (var J = 2; J <= N; J++) {
      dcopy(N - J + 1, AINV(JJ), 1, WORK(J, J - 1).asArray(), 1);
      dcopy(N - J, AINV(JJ + 1), 1, WORK(J, J).asArray(), LDWORK);
      JJ += N - J + 1;
    }

    // Multiply by A

    for (var J = N; J >= 2; J--) {
      dspmv('Lower', N, -ONE, A, WORK(1, J - 1).asArray(), 1, ZERO,
          WORK(1, J).asArray(), 1);
    }
    dspmv('Lower', N, -ONE, A, AINV(1), 1, ZERO, WORK(1, 1).asArray(), 1);
  }

  // Add the identity matrix to WORK .

  for (var I = 1; I <= N; I++) {
    WORK[I][I] += ONE;
  }

  // Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)

  RESID.value = dlange('1', N, N, WORK, LDWORK, RWORK);

  RESID.value = ((RESID.value * RCOND.value) / EPS) / N;
}
