// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dzasum.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlagtm.dart';
import 'package:dart_lapack/src/zlangt.dart';

void zgtt02(
  final String TRANS,
  final int N,
  final int NRHS,
  final Array<Complex> DL_,
  final Array<Complex> D_,
  final Array<Complex> DU_,
  final Matrix<Complex> X_,
  final int LDX,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<double> RESID,
) {
  final DL = DL_.having();
  final D = D_.having();
  final DU = DU_.having();
  final X = X_.having(ld: LDX);
  final B = B_.having(ld: LDB);
  const ONE = 1.0, ZERO = 0.0;

  // Quick exit if N = 0 or NRHS = 0

  RESID.value = ZERO;
  if (N <= 0 || NRHS == 0) return;

  // Compute the maximum over the number of right hand sides of
  //    norm(B - op(A)*X) / ( norm(op(A)) * norm(X) * EPS ).

  final double ANORM;
  if (lsame(TRANS, 'N')) {
    ANORM = zlangt('1', N, DL, D, DU);
  } else {
    ANORM = zlangt('I', N, DL, D, DU);
  }

  // Exit with RESID = 1/EPS if ANORM = 0.

  final EPS = dlamch('Epsilon');
  if (ANORM <= ZERO) {
    RESID.value = ONE / EPS;
    return;
  }

  // Compute B - op(A)*X and store in B.

  zlagtm(TRANS, N, NRHS, -ONE, DL, D, DU, X, LDX, ONE, B, LDB);

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
