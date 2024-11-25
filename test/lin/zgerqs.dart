// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';

void zgerqs(
  final int M,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> TAU_,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final TAU = TAU_.having();
  final WORK = WORK_.having(length: LWORK);

  // Test the input parameters.

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0 || M > N) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, M)) {
    INFO.value = -5;
  } else if (LDB < max(1, N)) {
    INFO.value = -8;
  } else if (LWORK < 1 || LWORK < NRHS && M > 0 && N > 0) {
    INFO.value = -10;
  }
  if (INFO.value != 0) {
    xerbla('ZGERQS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0 || M == 0) return;

  // Solve R*X = B(n-m+1:n,:)

  ztrsm('Left', 'Upper', 'No transpose', 'Non-unit', M, NRHS, Complex.one,
      A(1, N - M + 1), LDA, B(N - M + 1, 1), LDB);

  // Set B(1:n-m,:) to zero

  zlaset('Full', N - M, NRHS, Complex.zero, Complex.zero, B, LDB);

  // B := Q' * B

  zunmrq('Left', 'Conjugate transpose', N, NRHS, M, A, LDA, TAU, B, LDB, WORK,
      LWORK, INFO);
}
