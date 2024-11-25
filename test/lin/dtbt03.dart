// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/daxpy.dart';
import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/blas/dtbmv.dart';
import 'package:dart_lapack/src/blas/idamax.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';

void dtbt03(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int N,
  final int KD,
  final int NRHS,
  final Matrix<double> AB_,
  final int LDAB,
  final double SCALE,
  final Array<double> CNORM_,
  final double TSCAL,
  final Matrix<double> X_,
  final int LDX,
  final Matrix<double> B_,
  final int LDB,
  final Array<double> WORK_,
  final Box<double> RESID,
) {
  final AB = AB_.having(ld: LDAB);
  final X = X_.having(ld: LDX);
  final B = B_.having(ld: LDB);
  final CNORM = CNORM_.having();
  final WORK = WORK_.having();
  const ONE = 1.0, ZERO = 0.0;

  // Quick exit if N = 0

  if (N <= 0 || NRHS <= 0) {
    RESID.value = ZERO;
    return;
  }
  final EPS = dlamch('Epsilon');
  final SMLNUM = dlamch('Safe minimum');

  // Compute the norm of the triangular matrix A using the column
  // norms already computed by DLATBS.

  var TNORM = ZERO;
  if (lsame(DIAG, 'N')) {
    if (lsame(UPLO, 'U')) {
      for (var J = 1; J <= N; J++) {
        TNORM = max(TNORM, TSCAL * AB[KD + 1][J].abs() + CNORM[J]);
      }
    } else {
      for (var J = 1; J <= N; J++) {
        TNORM = max(TNORM, TSCAL * AB[1][J].abs() + CNORM[J]);
      }
    }
  } else {
    for (var J = 1; J <= N; J++) {
      TNORM = max(TNORM, TSCAL + CNORM[J]);
    }
  }

  // Compute the maximum over the number of right hand sides of
  //    norm(op(A)*x - s*b) / ( norm(op(A)) * norm(x) * EPS ).

  RESID.value = ZERO;
  for (var J = 1; J <= NRHS; J++) {
    dcopy(N, X(1, J).asArray(), 1, WORK, 1);
    var IX = idamax(N, WORK, 1);
    var XNORM = max(ONE, X[IX][J].abs());
    final XSCAL = (ONE / XNORM) / (KD + 1);
    dscal(N, XSCAL, WORK, 1);
    dtbmv(UPLO, TRANS, DIAG, N, KD, AB, LDAB, WORK, 1);
    daxpy(N, -SCALE * XSCAL, B(1, J).asArray(), 1, WORK, 1);
    IX = idamax(N, WORK, 1);
    var ERR = TSCAL * WORK[IX].abs();
    IX = idamax(N, X(1, J).asArray(), 1);
    XNORM = X[IX][J].abs();
    if (ERR * SMLNUM <= XNORM) {
      if (XNORM > ZERO) ERR /= XNORM;
    } else {
      if (ERR > ZERO) ERR = ONE / EPS;
    }
    if (ERR * SMLNUM <= TNORM) {
      if (TNORM > ZERO) ERR /= TNORM;
    } else {
      if (ERR > ZERO) ERR = ONE / EPS;
    }
    RESID.value = max(RESID.value, ERR);
  }
}
