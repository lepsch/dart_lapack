// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dzasum.dart';
import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlange.dart';

void zqrt16(
  final String TRANS,
  final int M,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> X_,
  final int LDX,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
  final A = A_.having(ld: LDA);
  final X = X_.having(ld: LDX);
  final B = B_.having(ld: LDB);
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if M = 0 or N = 0 or NRHS = 0

  if (M <= 0 || N <= 0 || NRHS == 0) {
    RESID.value = ZERO;
    return;
  }

  final double ANORM;
  final int N1, N2;
  if (lsame(TRANS, 'T') || lsame(TRANS, 'C')) {
    ANORM = zlange('I', M, N, A, LDA, RWORK);
    N1 = N;
    N2 = M;
  } else {
    ANORM = zlange('1', M, N, A, LDA, RWORK);
    N1 = M;
    N2 = N;
  }

  final EPS = dlamch('Epsilon');

  // Compute  B - A*X  (or  B - A'*X ) and store in B.

  zgemm(TRANS, 'No transpose', N1, NRHS, N2, -Complex.one, A, LDA, X, LDX,
      Complex.one, B, LDB);

  // Compute the maximum over the number of right hand sides of
  //    norm(B - A*X) / ( max(m,n) * norm(A) * norm(X) * EPS ) .

  RESID.value = ZERO;
  for (var J = 1; J <= NRHS; J++) {
    final BNORM = dzasum(N1, B(1, J).asArray(), 1);
    final XNORM = dzasum(N2, X(1, J).asArray(), 1);
    if (ANORM == ZERO && BNORM == ZERO) {
      RESID.value = ZERO;
    } else if (ANORM <= ZERO || XNORM <= ZERO) {
      RESID.value = ONE / EPS;
    } else {
      RESID.value =
          max(RESID.value, ((BNORM / ANORM) / XNORM) / (max(M, N) * EPS));
    }
  }
}
