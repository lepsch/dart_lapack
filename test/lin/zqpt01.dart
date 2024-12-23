// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zaxpy.dart';
import 'package:dart_lapack/src/blas/zcopy.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zunmqr.dart';

double zqpt01(
  final int M,
  final int N,
  final int K,
  final Matrix<Complex> A_,
  final Matrix<Complex> AF_,
  final int LDA,
  final Array<Complex> TAU_,
  final Array<int> JPVT_,
  final Array<Complex> WORK_,
  final int LWORK,
) {
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDA);
  final TAU = TAU_.having();
  final JPVT = JPVT_.having();
  final WORK = WORK_.having(length: LWORK);

  const ZERO = 0.0;
  final RWORK = Array<double>(1);

  // Test if there is enough workspace

  if (LWORK < M * N + N) {
    xerbla('ZQPT01', 10);
    return ZERO;
  }

  // Quick return if possible

  if (M <= 0 || N <= 0) return ZERO;

  final NORMA = zlange('One-norm', M, N, A, LDA, RWORK);

  for (var J = 1; J <= K; J++) {
    for (var I = 1; I <= min(J, M); I++) {
      WORK[(J - 1) * M + I] = AF[I][J];
    }
    for (var I = J + 1; I <= M; I++) {
      WORK[(J - 1) * M + I] = Complex.zero;
    }
  }
  for (var J = K + 1; J <= N; J++) {
    zcopy(M, AF(1, J).asArray(), 1, WORK((J - 1) * M + 1), 1);
  }

  zunmqr('Left', 'No transpose', M, N, K, AF, LDA, TAU, WORK.asMatrix(), M,
      WORK(M * N + 1), LWORK - M * N, Box(0));

  for (var J = 1; J <= N; J++) {
    // Compare i-th column of QR and jpvt(i)-th column of A

    zaxpy(
        M, -Complex.one, A(1, JPVT[J]).asArray(), 1, WORK((J - 1) * M + 1), 1);
  }

  final result = zlange('One-norm', M, N, WORK.asMatrix(), M, RWORK) /
      (max(M, N) * dlamch('Epsilon'));
  return NORMA != ZERO ? result / NORMA : result;
}
