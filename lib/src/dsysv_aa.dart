// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dsytrf_aa.dart';
import 'package:dart_lapack/src/dsytrs_aa.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dsysv_aa(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Matrix<double> B_,
  final int LDB,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  final WORK = WORK_.having();
  bool LQUERY;
  int LWKMIN, LWKOPT = 0, LWKOPT_SYTRF, LWKOPT_SYTRS;

  // Test the input parameters.

  INFO.value = 0;
  LQUERY = (LWORK == -1);
  LWKMIN = max(1, max(2 * N, 3 * N - 2));
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDB < max(1, N)) {
    INFO.value = -8;
  } else if (LWORK < LWKMIN && !LQUERY) {
    INFO.value = -10;
  }

  if (INFO.value == 0) {
    dsytrf_aa(UPLO, N, A, LDA, IPIV, WORK, -1, INFO);
    LWKOPT_SYTRF = WORK[1].toInt();
    dsytrs_aa(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, -1, INFO);
    LWKOPT_SYTRS = WORK[1].toInt();
    LWKOPT = max(LWKMIN, max(LWKOPT_SYTRF, LWKOPT_SYTRS));
    WORK[1] = LWKOPT.toDouble();
  }

  if (INFO.value != 0) {
    xerbla('DSYSV_AA', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Compute the factorization A = U**T*T*U or A = L*T*L**T.

  dsytrf_aa(UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO);
  if (INFO.value == 0) {
    // Solve the system A*X = B, overwriting B with X.

    dsytrs_aa(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO);
  }

  WORK[1] = LWKOPT.toDouble();
}
