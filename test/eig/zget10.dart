// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dzasum.dart';
import 'package:dart_lapack/src/blas/zaxpy.dart';
import 'package:dart_lapack/src/blas/zcopy.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlange.dart';

void zget10(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Box<double> RESULT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  int J;
  double ANORM, EPS, UNFL, WNORM;

  // Quick return if possible

  if (M <= 0 || N <= 0) {
    RESULT.value = ZERO;
    return;
  }

  UNFL = dlamch('Safe minimum');
  EPS = dlamch('Precision');

  WNORM = ZERO;
  for (J = 1; J <= N; J++) {
    zcopy(M, A(1, J).asArray(), 1, WORK, 1);
    zaxpy(M, -Complex.one, B(1, J).asArray(), 1, WORK, 1);
    WNORM = max(WNORM, dzasum(N, WORK, 1));
  }

  ANORM = max(zlange('1', M, N, A, LDA, RWORK), UNFL);

  if (ANORM > WNORM) {
    RESULT.value = (WNORM / ANORM) / (M * EPS);
  } else {
    if (ANORM < ONE) {
      RESULT.value = (min(WNORM, M * ANORM) / ANORM) / (M * EPS);
    } else {
      RESULT.value = min(WNORM / ANORM, M) / (M * EPS);
    }
  }
}
