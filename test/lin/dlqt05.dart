// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';

void dlqt05(
  final int M,
  final int N,
  final int L,
  final int NB,
  final Array<double> RESULT_,
) {
  final RESULT = RESULT_.having(length: 6);
  const ZERO = 0.0, ONE = 1.0;
  final ISEED = Array.fromList([1988, 1989, 1990, 1991]);
  final INFO = Box(0);

  final EPS = dlamch('Epsilon');
  final K = M;
  final N2 = M + N;
  final NP1 = N > 0 ? M + 1 : 1;
  final LWORK = N2 * N2 * NB;

  final A = Matrix<double>(M, N2),
      AF = Matrix<double>(M, N2),
      Q = Matrix<double>(N2, N2),
      R = Matrix<double>(N2, N2),
      RWORK = Array<double>(N2),
      WORK = Array<double>(LWORK),
      T = Matrix<double>(NB, M),
      C = Matrix<double>(N2, M),
      CF = Matrix<double>(N2, M),
      D = Matrix<double>(M, N2),
      DF = Matrix<double>(M, N2);

  // Put random stuff into A

  final LDT = NB;
  dlaset('Full', M, N2, ZERO, ZERO, A, M);
  dlaset('Full', NB, M, ZERO, ZERO, T, NB);
  for (var J = 1; J <= M; J++) {
    dlarnv(2, ISEED, M - J + 1, A(J, J).asArray());
  }
  if (N > 0) {
    for (var J = 1; J <= N - L; J++) {
      dlarnv(2, ISEED, M, A(1, min(N + M, M + 1) + J - 1).asArray());
    }
  }
  if (L > 0) {
    for (var J = 1; J <= L; J++) {
      dlarnv(2, ISEED, M - J + 1,
          A(J, min(N + M, N + M - L + 1) + J - 1).asArray());
    }
  }

  // Copy the matrix A to the array AF.

  dlacpy('Full', M, N2, A, M, AF, M);

  // Factor the matrix A in the array AF.

  dtplqt(M, N, L, NB, AF, M, AF(1, NP1), M, T, LDT, WORK, INFO);

  // Generate the (M+N)-by-(M+N) matrix Q by applying H to I

  dlaset('Full', N2, N2, ZERO, ONE, Q, N2);
  dgemlqt('L', 'N', N2, N2, K, NB, AF, M, T, LDT, Q, N2, WORK, INFO);

  // Copy L

  dlaset('Full', N2, N2, ZERO, ZERO, R, N2);
  dlacpy('Lower', M, N2, AF, M, R, N2);

  // Compute |L - A*Q*T| / |A| and store in RESULT(1)

  dgemm('N', 'T', M, N2, N2, -ONE, A, M, Q, N2, ONE, R, N2);
  final ANORM = dlange('1', M, N2, A, M, RWORK);
  var RESID = dlange('1', M, N2, R, N2, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = RESID / (EPS * ANORM * max(1, N2));
  } else {
    RESULT[1] = ZERO;
  }

  // Compute |I - Q*Q'| and store in RESULT(2)

  dlaset('Full', N2, N2, ZERO, ONE, R, N2);
  dsyrk('U', 'N', N2, N2, -ONE, Q, N2, ONE, R, N2);
  RESID = dlansy('1', 'Upper', N2, R, N2, RWORK);
  RESULT[2] = RESID / (EPS * max(1, N2));

  // Generate random m-by-n matrix C and a copy CF

  dlaset('Full', N2, M, ZERO, ONE, C, N2);
  for (var J = 1; J <= M; J++) {
    dlarnv(2, ISEED, N2, C(1, J).asArray());
  }
  final CNORM = dlange('1', N2, M, C, N2, RWORK);
  dlacpy('Full', N2, M, C, N2, CF, N2);

  // Apply Q to C as Q*C

  dtpmlqt('L', 'N', N, M, K, L, NB, AF(1, NP1), M, T, LDT, CF, N2, CF(NP1, 1),
      N2, WORK, INFO);

  // Compute |Q*C - Q*C| / |C|

  dgemm('N', 'N', N2, M, N2, -ONE, Q, N2, C, N2, ONE, CF, N2);
  RESID = dlange('1', N2, M, CF, N2, RWORK);
  if (CNORM > ZERO) {
    RESULT[3] = RESID / (EPS * max(1, N2) * CNORM);
  } else {
    RESULT[3] = ZERO;
  }

  // Copy C into CF again

  dlacpy('Full', N2, M, C, N2, CF, N2);

  // Apply Q to C as QT*C

  dtpmlqt('L', 'T', N, M, K, L, NB, AF(1, NP1), M, T, LDT, CF, N2, CF(NP1, 1),
      N2, WORK, INFO);

  // Compute |QT*C - QT*C| / |C|

  dgemm('T', 'N', N2, M, N2, -ONE, Q, N2, C, N2, ONE, CF, N2);
  RESID = dlange('1', N2, M, CF, N2, RWORK);

  if (CNORM > ZERO) {
    RESULT[4] = RESID / (EPS * max(1, N2) * CNORM);
  } else {
    RESULT[4] = ZERO;
  }

  // Generate random m-by-n matrix D and a copy DF

  for (var J = 1; J <= N2; J++) {
    dlarnv(2, ISEED, M, D(1, J).asArray());
  }
  final DNORM = dlange('1', M, N2, D, M, RWORK);
  dlacpy('Full', M, N2, D, M, DF, M);

  // Apply Q to D as D*Q

  dtpmlqt('R', 'N', M, N, K, L, NB, AF(1, NP1), M, T, LDT, DF, M, DF(1, NP1), M,
      WORK, INFO);

  // Compute |D*Q - D*Q| / |D|

  dgemm('N', 'N', M, N2, N2, -ONE, D, M, Q, N2, ONE, DF, M);
  RESID = dlange('1', M, N2, DF, M, RWORK);
  if (CNORM > ZERO) {
    RESULT[5] = RESID / (EPS * max(1, N2) * DNORM);
  } else {
    RESULT[5] = ZERO;
  }

  // Copy D into DF again

  dlacpy('Full', M, N2, D, M, DF, M);

  // Apply Q to D as D*QT

  dtpmlqt('R', 'T', M, N, K, L, NB, AF(1, NP1), M, T, LDT, DF, M, DF(1, NP1), M,
      WORK, INFO);

  // Compute |D*QT - D*QT| / |D|

  dgemm('N', 'T', M, N2, N2, -ONE, D, M, Q, N2, ONE, DF, M);
  RESID = dlange('1', M, N2, DF, M, RWORK);
  if (CNORM > ZERO) {
    RESULT[6] = RESID / (EPS * max(1, N2) * DNORM);
  } else {
    RESULT[6] = ZERO;
  }
}
