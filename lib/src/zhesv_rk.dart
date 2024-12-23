// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zhetrf_rk.dart';
import 'package:dart_lapack/src/zhetrs_3.dart';

void zhesv_rk(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> E_,
  final Array<int> IPIV_,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  final E = E_.having();
  final WORK = WORK_.having();
  bool LQUERY;
  int LWKOPT = 0;

  // Test the input parameters.

  INFO.value = 0;
  LQUERY = (LWORK == -1);
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDB < max(1, N)) {
    INFO.value = -9;
  } else if (LWORK < 1 && !LQUERY) {
    INFO.value = -11;
  }

  if (INFO.value == 0) {
    if (N == 0) {
      LWKOPT = 1;
    } else {
      zhetrf_rk(UPLO, N, A, LDA, E, IPIV, WORK, -1, INFO);
      LWKOPT = WORK[1].toInt();
    }
    WORK[1] = LWKOPT.toComplex();
  }

  if (INFO.value != 0) {
    xerbla('ZHESV_RK', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Compute the factorization A = P*U*D*(U**H)*(P**T) or
  // A = P*U*D*(U**H)*(P**T).

  zhetrf_rk(UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO);

  if (INFO.value == 0) {
    // Solve the system A*X = B with BLAS3 solver, overwriting B with X.

    zhetrs_3(UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, INFO);
  }

  WORK[1] = LWKOPT.toComplex();
}
