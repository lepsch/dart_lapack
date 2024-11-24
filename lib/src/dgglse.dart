// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dtrmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dggrqf.dart';
import 'package:lapack/src/dormqr.dart';
import 'package:lapack/src/dormrq.dart';
import 'package:lapack/src/dtrtrs.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgglse(
  final int M,
  final int N,
  final int P,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Array<double> C_,
  final Array<double> D_,
  final Array<double> X_,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final C = C_.having();
  final D = D_.having();
  final X = X_.having();
  final WORK = WORK_.having();
  const ONE = 1.0;
  bool LQUERY;
  int LOPT, LWKMIN, LWKOPT, MN, NB, NB1, NB2, NB3, NB4, NR;

  // Test the input parameters

  INFO.value = 0;
  MN = min(M, N);
  LQUERY = (LWORK == -1);
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (P < 0 || P > N || P < N - M) {
    INFO.value = -3;
  } else if (LDA < max(1, M)) {
    INFO.value = -5;
  } else if (LDB < max(1, P)) {
    INFO.value = -7;
  }

  // Calculate workspace

  if (INFO.value == 0) {
    if (N == 0) {
      LWKMIN = 1;
      LWKOPT = 1;
    } else {
      NB1 = ilaenv(1, 'DGEQRF', ' ', M, N, -1, -1);
      NB2 = ilaenv(1, 'DGERQF', ' ', M, N, -1, -1);
      NB3 = ilaenv(1, 'DORMQR', ' ', M, N, P, -1);
      NB4 = ilaenv(1, 'DORMRQ', ' ', M, N, P, -1);
      NB = max(max(NB1, NB2), max(NB3, NB4));
      LWKMIN = M + N + P;
      LWKOPT = P + MN + max(M, N) * NB;
    }
    WORK[1] = LWKOPT.toDouble();

    if (LWORK < LWKMIN && !LQUERY) {
      INFO.value = -12;
    }
  }

  if (INFO.value != 0) {
    xerbla('DGGLSE', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Compute the GRQ factorization of matrices B and A:

  // B*Q**T = (  0  T12 ) P   Z**T*A*Q**T = ( R11 R12 ) N-P
  //             N-P  P                     (  0  R22 ) M+P-N
  //                                           N-P  P

  // where T12 and R11 are upper triangular, and Q and Z are
  // orthogonal.

  dggrqf(P, M, N, B, LDB, WORK, A, LDA, WORK(P + 1), WORK(P + MN + 1),
      LWORK - P - MN, INFO);
  LOPT = WORK[P + MN + 1].toInt();

  // Update c = Z**T *c = ( c1 ) N-P
  //                      ( c2 ) M+P-N

  dormqr('Left', 'Transpose', M, 1, MN, A, LDA, WORK(P + 1),
      C.asMatrix(max(1, M)), max(1, M), WORK(P + MN + 1), LWORK - P - MN, INFO);
  LOPT = max(LOPT, WORK[P + MN + 1].toInt());

  // Solve T12*x2 = d for x2

  if (P > 0) {
    dtrtrs('Upper', 'No transpose', 'Non-unit', P, 1, B(1, N - P + 1), LDB,
        D.asMatrix(P), P, INFO);

    if (INFO.value > 0) {
      INFO.value = 1;
      return;
    }

    // Put the solution in X

    dcopy(P, D, 1, X(N - P + 1), 1);

    // Update c1

    dgemv(
        'No transpose', N - P, P, -ONE, A(1, N - P + 1), LDA, D, 1, ONE, C, 1);
  }

  // Solve R11*x1 = c1 for x1

  if (N > P) {
    dtrtrs('Upper', 'No transpose', 'Non-unit', N - P, 1, A, LDA,
        C.asMatrix(N - P), N - P, INFO);

    if (INFO.value > 0) {
      INFO.value = 2;
      return;
    }

    // Put the solutions in X

    dcopy(N - P, C, 1, X, 1);
  }

  // Compute the residual vector:

  if (M < N) {
    NR = M + P - N;
    if (NR > 0) {
      dgemv('No transpose', NR, N - M, -ONE, A(N - P + 1, M + 1), LDA,
          D(NR + 1), 1, ONE, C(N - P + 1), 1);
    }
  } else {
    NR = P;
  }
  if (NR > 0) {
    dtrmv('Upper', 'No transpose', 'Non unit', NR, A(N - P + 1, N - P + 1), LDA,
        D, 1);
    daxpy(NR, -ONE, D, 1, C(N - P + 1), 1);
  }

  // Backward transformation x = Q**T*x

  dormrq('Left', 'Transpose', N, 1, P, B, LDB, WORK(1), X.asMatrix(N), N,
      WORK(P + MN + 1), LWORK - P - MN, INFO);
  WORK[1] = P + MN + max(LOPT, WORK[P + MN + 1].toInt()).toDouble();
}
