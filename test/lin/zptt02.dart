// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dzasum.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlanht.dart';

import 'zlaptm.dart';

void zptt02(
  final String UPLO,
  final int N,
  final int NRHS,
  final Array<double> D_,
  final Array<Complex> E_,
  final Matrix<Complex> X_,
  final int LDX,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<double> RESID,
) {
  final D = D_.having();
  final E = E_.having();
  final X = X_.having(ld: LDX);
  final B = B_.having(ld: LDB);
  const ONE = 1.0, ZERO = 0.0;

  // Quick return if possible

  if (N <= 0) {
    RESID.value = ZERO;
    return;
  }

  // Compute the 1-norm of the tridiagonal matrix A.

  final ANORM = zlanht('1', N, D, E);

  // Exit with RESID = 1/EPS if ANORM = 0.

  final EPS = dlamch('Epsilon');
  if (ANORM <= ZERO) {
    RESID.value = ONE / EPS;
    return;
  }

  // Compute B - A*X.

  zlaptm(UPLO, N, NRHS, -ONE, D, E, X, LDX, ONE, B, LDB);

  // Compute the maximum over the number of right hand sides of
  //    norm(B - A*X) / ( norm(A) * norm(X) * EPS ).

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
