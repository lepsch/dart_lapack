// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dzasum.dart';
import 'package:dart_lapack/src/blas/zhpmv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlanhp.dart';

void zppt02(
  final String UPLO,
  final int N,
  final int NRHS,
  final Array<Complex> A_,
  final Matrix<Complex> X_,
  final int LDX,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
  final A = A_.having();
  final X = X_.having(ld: LDX);
  final B = B_.having(ld: LDB);
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if N = 0 or NRHS = 0.

  if (N <= 0 || NRHS <= 0) {
    RESID.value = ZERO;
    return;
  }

  // Exit with RESID = 1/EPS if ANORM = 0.

  final EPS = dlamch('Epsilon');
  final ANORM = zlanhp('1', UPLO, N, A, RWORK);
  if (ANORM <= ZERO) {
    RESID.value = ONE / EPS;
    return;
  }

  // Compute  B - A*X  for the matrix of right hand sides B.

  for (var J = 1; J <= NRHS; J++) {
    zhpmv(UPLO, N, -Complex.one, A, X(1, J).asArray(), 1, Complex.one,
        B(1, J).asArray(), 1);
  }

  // Compute the maximum over the number of right hand sides of
  //    norm( B - A*X ) / ( norm(A) * norm(X) * EPS ) .

  RESID.value = ZERO;
  for (var J = 1; J <= NRHS; J++) {
    final BNORM = dzasum(N, B(1, J).asArray(), 1);
    final XNORM = dzasum(N, X(1, J).asArray(), 1);
    if (XNORM <= ZERO) {
      RESID.value = ONE / EPS;
    } else {
      RESID.value = max(RESID.value, ((BNORM / ANORM) / XNORM) / EPS);
    }
  }
}
