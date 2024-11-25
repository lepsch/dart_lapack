// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dzasum.dart';
import 'package:dart_lapack/src/blas/zaxpy.dart';
import 'package:dart_lapack/src/blas/zcopy.dart';
import 'package:dart_lapack/src/blas/ztpmv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlantp.dart';

void ztpt02(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int N,
  final int NRHS,
  final Array<Complex> AP_,
  final Matrix<Complex> X_,
  final int LDX,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
  final AP = AP_.having();
  final X = X_.having(ld: LDX);
  final B = B_.having(ld: LDB);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if N = 0 or NRHS = 0

  if (N <= 0 || NRHS <= 0) {
    RESID.value = ZERO;
    return;
  }

  // Compute the 1-norm of op(A).

  final double ANORM;
  if (lsame(TRANS, 'N')) {
    ANORM = zlantp('1', UPLO, DIAG, N, AP, RWORK);
  } else {
    ANORM = zlantp('I', UPLO, DIAG, N, AP, RWORK);
  }

  // Exit with RESID = 1/EPS if ANORM = 0.

  final EPS = dlamch('Epsilon');
  if (ANORM <= ZERO) {
    RESID.value = ONE / EPS;
    return;
  }

  // Compute the maximum over the number of right hand sides of
  //    norm(op(A)*X - B) / ( norm(op(A)) * norm(X) * EPS ).

  RESID.value = ZERO;
  for (var J = 1; J <= NRHS; J++) {
    zcopy(N, X(1, J).asArray(), 1, WORK, 1);
    ztpmv(UPLO, TRANS, DIAG, N, AP, WORK, 1);
    zaxpy(N, Complex(-ONE), B(1, J).asArray(), 1, WORK, 1);
    final BNORM = dzasum(N, WORK, 1);
    final XNORM = dzasum(N, X(1, J).asArray(), 1);
    if (XNORM <= ZERO) {
      RESID.value = ONE / EPS;
    } else {
      RESID.value = max(RESID.value, ((BNORM / ANORM) / XNORM) / EPS);
    }
  }
}
