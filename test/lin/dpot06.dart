// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dsymm.dart';
import 'package:dart_lapack/src/blas/idamax.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlansy.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void dpot06(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> X_,
  final int LDX,
  final Matrix<double> B_,
  final int LDB,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
  final A = A_.having(ld: LDA);
  final X = X_.having(ld: LDX);
  final B = B_.having(ld: LDB);
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if N = 0 or NRHS = 0

  if (N <= 0 || NRHS == 0) {
    RESID.value = ZERO;
    return;
  }

  // Exit with RESID = 1/EPS if ANORM = 0.

  final EPS = dlamch('Epsilon');
  final ANORM = dlansy('I', UPLO, N, A, LDA, RWORK);
  if (ANORM <= ZERO) {
    RESID.value = ONE / EPS;
    return;
  }

  // Compute  B - A*X  and store in B.

  dsymm('Left', UPLO, N, NRHS, -ONE, A, LDA, X, LDX, ONE, B, LDB);

  // Compute the maximum over the number of right hand sides of
  //    norm(B - A*X) / ( norm(A) * norm(X) * EPS ) .

  RESID.value = ZERO;
  for (var J = 1; J <= NRHS; J++) {
    final BNORM = B[idamax(N, B(1, J).asArray(), 1)][J].abs();
    final XNORM = X[idamax(N, X(1, J).asArray(), 1)][J].abs();
    if (XNORM <= ZERO) {
      RESID.value = ONE / EPS;
    } else {
      RESID.value = max(RESID.value, ((BNORM / ANORM) / XNORM) / EPS);
    }
  }
}
