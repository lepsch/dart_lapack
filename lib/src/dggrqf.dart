// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgerqf.dart';
import 'package:dart_lapack/src/dormrq.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/dgeqrf.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dggrqf(
  final int M,
  final int P,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> TAUA_,
  final Matrix<double> B_,
  final int LDB,
  final Array<double> TAUB_,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final TAUA = TAUA_.having();
  final B = B_.having(ld: LDB);
  final TAUB = TAUB_.having();
  final WORK = WORK_.having();
  bool LQUERY;
  int LOPT, LWKOPT, NB, NB1, NB2, NB3;

  // Test the input parameters

  INFO.value = 0;
  NB1 = ilaenv(1, 'DGERQF', ' ', M, N, -1, -1);
  NB2 = ilaenv(1, 'DGEQRF', ' ', P, N, -1, -1);
  NB3 = ilaenv(1, 'DORMRQ', ' ', M, N, P, -1);
  NB = max(NB1, max(NB2, NB3));
  LWKOPT = max(1, max(N, max(M, P)) * NB);
  WORK[1] = LWKOPT.toDouble();
  LQUERY = (LWORK == -1);
  if (M < 0) {
    INFO.value = -1;
  } else if (P < 0) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, M)) {
    INFO.value = -5;
  } else if (LDB < max(1, P)) {
    INFO.value = -8;
  } else if (LWORK < max(max(1, M), max(P, N)) && !LQUERY) {
    INFO.value = -11;
  }
  if (INFO.value != 0) {
    xerbla('DGGRQF', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // RQ factorization of M-by-N matrix A: A = R*Q

  dgerqf(M, N, A, LDA, TAUA, WORK, LWORK, INFO);
  LOPT = WORK[1].toInt();

  // Update B := B*Q**T

  dormrq('Right', 'Transpose', P, N, min(M, N), A(max(1, M - N + 1), 1), LDA,
      TAUA, B, LDB, WORK, LWORK, INFO);
  LOPT = max(LOPT, WORK[1].toInt());

  // QR factorization of P-by-N matrix B: B = Z*T

  dgeqrf(P, N, B, LDB, TAUB, WORK, LWORK, INFO);
  WORK[1] = max(LOPT, WORK[1].toInt()).toDouble();
}
