// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlaset.dart';
import 'package:dart_lapack/src/zunm2r.dart';

double zqrt11(
  final int M,
  final int K,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final int LWORK,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having(length: LWORK);
  final INFO = Box(0);

  // Test for sufficient workspace

  if (LWORK < M * M + M) {
    xerbla('ZQRT11', 7);
    return 0;
  }

  // Quick return if possible

  if (M <= 0) return 0;

  zlaset('Full', M, M, Complex.zero, Complex.one, WORK.asMatrix(), M);

  // Form Q

  zunm2r('Left', 'No transpose', M, M, K, A, LDA, TAU, WORK.asMatrix(), M,
      WORK(M * M + 1), INFO);

  // Form Q'*Q

  zunm2r('Left', 'Conjugate transpose', M, M, K, A, LDA, TAU, WORK.asMatrix(),
      M, WORK(M * M + 1), INFO);

  for (var J = 1; J <= M; J++) {
    WORK[(J - 1) * M + J] -= Complex.one;
  }

  final RDUMMY = Array<double>(1);
  return zlange('One-norm', M, M, WORK.asMatrix(), M, RDUMMY) /
      (M * dlamch('Epsilon'));
}
