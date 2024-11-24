// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';

double zrzt02(
  final int M,
  final int N,
  final Matrix<Complex> AF_,
  final int LDA,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final int LWORK,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AF = AF_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having(length: LWORK);
  const ZERO = 0.0;
  final INFO = Box(0);

  if (LWORK < N * N + N) {
    xerbla('ZRZT02', 7);
    return ZERO;
  }

  // Quick return if possible

  if (M <= 0 || N <= 0) return ZERO;

  // Q := I

  zlaset('Full', N, N, Complex.zero, Complex.one, WORK.asMatrix(), N);

  // Q := P(1) * ... * P(m) * Q

  zunmrz('Left', 'No transpose', N, N, M, N - M, AF, LDA, TAU, WORK.asMatrix(),
      N, WORK(N * N + 1), LWORK - N * N, INFO);

  // Q := P(m)' * ... * P(1)' * Q

  zunmrz('Left', 'Conjugate transpose', N, N, M, N - M, AF, LDA, TAU,
      WORK.asMatrix(), N, WORK(N * N + 1), LWORK - N * N, INFO);

  // Q := Q - I

  for (var I = 1; I <= N; I++) {
    WORK[(I - 1) * N + I] -= Complex.one;
  }

  final RWORK = Array<double>(1);
  return zlange('One-norm', N, N, WORK.asMatrix(), N, RWORK) /
      (dlamch('Epsilon') * max(M, N));
}
