// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/izamax.dart';
import 'package:dart_lapack/src/blas/zaxpy.dart';
import 'package:dart_lapack/src/blas/zcopy.dart';
import 'package:dart_lapack/src/blas/zdscal.dart';
import 'package:dart_lapack/src/blas/ztrmv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';

void ztrt03(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final double SCALE,
  final Array<double> CNORM_,
  final double TSCAL,
  final Matrix<Complex> X_,
  final int LDX,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<Complex> WORK_,
  final Box<double> RESID,
) {
  final A = A_.having(ld: LDA);
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
  // norms already computed by ZLATRS.

  var TNORM = ZERO;
  if (lsame(DIAG, 'N')) {
    for (var J = 1; J <= N; J++) {
      TNORM = max(TNORM, TSCAL * A[J][J].abs() + CNORM[J]);
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
    zcopy(N, X(1, J).asArray(), 1, WORK, 1);
    var IX = izamax(N, WORK, 1);
    var XNORM = max(ONE, X[IX][J].abs());
    final XSCAL = (ONE / XNORM) / N;
    zdscal(N, XSCAL, WORK, 1);
    ztrmv(UPLO, TRANS, DIAG, N, A, LDA, WORK, 1);
    zaxpy(N, Complex(-SCALE * XSCAL), B(1, J).asArray(), 1, WORK, 1);
    IX = izamax(N, WORK, 1);
    var ERR = TSCAL * WORK[IX].abs();
    IX = izamax(N, X(1, J).asArray(), 1);
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
