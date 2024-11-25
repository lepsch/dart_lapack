// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dzasum.dart';
import 'package:dart_lapack/src/blas/zcopy.dart';
import 'package:dart_lapack/src/blas/zgemv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zggglm.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';

void zglmts(
  final int N,
  final int M,
  final int P,
  final Matrix<Complex> A_,
  final Matrix<Complex> AF_,
  final int LDA,
  final Matrix<Complex> B_,
  final Matrix<Complex> BF_,
  final int LDB,
  final Array<Complex> D_,
  final Array<Complex> DF_,
  final Array<Complex> X_,
  final Array<Complex> U_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Box<double> RESULT,
) {
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final BF = BF_.having(ld: LDB);
  final D = D_.having();
  final DF = DF_.having();
  final X = X_.having();
  final U = U_.having();
  final WORK = WORK_.having(length: LWORK);
  final RWORK = RWORK_.having();
  const ZERO = 0.0;
  final INFO = Box(0);
  double ANORM, BNORM, DNORM, EPS, UNFL, XNORM, YNORM;

  EPS = dlamch('Epsilon');
  UNFL = dlamch('Safe minimum');
  ANORM = max(zlange('1', N, M, A, LDA, RWORK), UNFL);
  BNORM = max(zlange('1', N, P, B, LDB, RWORK), UNFL);

  // Copy the matrices A and B to the arrays AF and BF,
  // and the vector D the array DF.

  zlacpy('Full', N, M, A, LDA, AF, LDA);
  zlacpy('Full', N, P, B, LDB, BF, LDB);
  zcopy(N, D, 1, DF, 1);

  // Solve GLM problem

  zggglm(N, M, P, AF, LDA, BF, LDB, DF, X, U, WORK, LWORK, INFO);

  // Test the residual for the solution of LSE

  //                   norm( d - A*x - B*u )
  //   RESULT = -----------------------------------------
  //            (norm(A)+norm(B))*(norm(x)+norm(u))*EPS

  zcopy(N, D, 1, DF, 1);
  zgemv('No transpose', N, M, -Complex.one, A, LDA, X, 1, Complex.one, DF, 1);

  zgemv('No transpose', N, P, -Complex.one, B, LDB, U, 1, Complex.one, DF, 1);

  DNORM = dzasum(N, DF, 1);
  XNORM = dzasum(M, X, 1) + dzasum(P, U, 1);
  YNORM = ANORM + BNORM;

  if (XNORM <= ZERO) {
    RESULT.value = ZERO;
  } else {
    RESULT.value = ((DNORM / YNORM) / XNORM) / EPS;
  }
}
