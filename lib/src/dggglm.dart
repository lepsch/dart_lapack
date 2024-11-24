// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/dgemv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dggqrf.dart';
import 'package:dart_lapack/src/dormqr.dart';
import 'package:dart_lapack/src/dormrq.dart';
import 'package:dart_lapack/src/dtrtrs.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dggglm(
  final int N,
  final int M,
  final int P,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Array<double> D_,
  final Array<double> X_,
  final Array<double> Y_,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final D = D_.having();
  final X = X_.having();
  final Y = Y_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool LQUERY;
  int I, LOPT = 0, LWKMIN, LWKOPT, NB, NB1, NB2, NB3, NB4, NP;

  // Test the input parameters

  INFO.value = 0;
  NP = min(N, P);
  LQUERY = (LWORK == -1);
  if (N < 0) {
    INFO.value = -1;
  } else if (M < 0 || M > N) {
    INFO.value = -2;
  } else if (P < 0 || P < N - M) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDB < max(1, N)) {
    INFO.value = -7;
  }

  // Calculate workspace

  if (INFO.value == 0) {
    if (N == 0) {
      LWKMIN = 1;
      LWKOPT = 1;
    } else {
      NB1 = ilaenv(1, 'DGEQRF', ' ', N, M, -1, -1);
      NB2 = ilaenv(1, 'DGERQF', ' ', N, M, -1, -1);
      NB3 = ilaenv(1, 'DORMQR', ' ', N, M, P, -1);
      NB4 = ilaenv(1, 'DORMRQ', ' ', N, M, P, -1);
      NB = max(max(NB1, NB2), max(NB3, NB4));
      LWKMIN = M + N + P;
      LWKOPT = M + NP + max(N, P) * NB;
    }
    WORK[1] = LWKOPT.toDouble();

    if (LWORK < LWKMIN && !LQUERY) {
      INFO.value = -12;
    }
  }

  if (INFO.value != 0) {
    xerbla('DGGGLM', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) {
    for (I = 1; I <= M; I++) {
      X[I] = ZERO;
    }
    for (I = 1; I <= P; I++) {
      Y[I] = ZERO;
    }
    return;
  }

  // Compute the GQR factorization of matrices A and B:

  // Q**T*A = ( R11 ) M,    Q**T*B*Z**T = ( T11   T12 ) M
  //          (  0  ) N-M                 (  0    T22 ) N-M
  //             M                         M+P-N  N-M

  // where R11 and T22 are upper triangular, and Q and Z are
  // orthogonal.

  dggqrf(N, M, P, A, LDA, WORK, B, LDB, WORK(M + 1), WORK(M + NP + 1),
      LWORK - M - NP, INFO);
  LOPT = WORK[M + NP + 1].toInt();

  // Update left-hand-side vector d = Q**T*d = ( d1 ) M
  //                                           ( d2 ) N-M

  dormqr('Left', 'Transpose', N, 1, M, A, LDA, WORK, D.asMatrix(max(1, N)),
      max(1, N), WORK(M + NP + 1), LWORK - M - NP, INFO);
  LOPT = max(LOPT, WORK[M + NP + 1].toInt());

  // Solve T22*y2 = d2 for y2

  if (N > M) {
    dtrtrs('Upper', 'No transpose', 'Non unit', N - M, 1,
        B(M + 1, M + P - N + 1), LDB, D(M + 1).asMatrix(N - M), N - M, INFO);

    if (INFO.value > 0) {
      INFO.value = 1;
      return;
    }

    dcopy(N - M, D(M + 1), 1, Y(M + P - N + 1), 1);
  }

  // Set y1 = 0

  for (I = 1; I <= M + P - N; I++) {
    Y[I] = ZERO;
  }

  // Update d1 -= T12*y2

  dgemv('No transpose', M, N - M, -ONE, B(1, M + P - N + 1), LDB,
      Y(M + P - N + 1), 1, ONE, D, 1);

  // Solve triangular system: R11*x = d1

  if (M > 0) {
    dtrtrs('Upper', 'No Transpose', 'Non unit', M, 1, A, LDA, D.asMatrix(M), M,
        INFO);

    if (INFO.value > 0) {
      INFO.value = 2;
      return;
    }

    // Copy D to X

    dcopy(M, D, 1, X, 1);
  }

  // Backward transformation y = Z**T *y

  dormrq(
      'Left',
      'Transpose',
      P,
      1,
      NP,
      B(max(1, N - P + 1), 1),
      LDB,
      WORK(M + 1),
      Y.asMatrix(max(1, P)),
      max(1, P),
      WORK(M + NP + 1),
      LWORK - M - NP,
      INFO);
  WORK[1] = M + NP + max(LOPT, WORK[M + NP + 1].toInt()).toDouble();
}
