// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zhetrf.dart';
import 'package:dart_lapack/src/zhetrs.dart';
import 'package:dart_lapack/src/zhetrs2.dart';

void zhesv(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
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
  final WORK = WORK_.having();
  bool LQUERY;
  int LWKOPT = 0, NB;

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
    INFO.value = -8;
  } else if (LWORK < 1 && !LQUERY) {
    INFO.value = -10;
  }

  if (INFO.value == 0) {
    if (N == 0) {
      LWKOPT = 1;
    } else {
      NB = ilaenv(1, 'ZHETRF', UPLO, N, -1, -1, -1);
      LWKOPT = N * NB;
    }
    WORK[1] = LWKOPT.toComplex();
  }

  if (INFO.value != 0) {
    xerbla('ZHESV', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Compute the factorization A = U*D*U**H or A = L*D*L**H.

  zhetrf(UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO);
  if (INFO.value == 0) {
    // Solve the system A*X = B, overwriting B with X.

    if (LWORK < N) {
      // Solve with TRS ( Use Level BLAS 2)

      zhetrs(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO);
    } else {
      // Solve with TRS2 ( Use Level BLAS 3)

      zhetrs2(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, INFO);
    }
  }

  WORK[1] = LWKOPT.toComplex();
}
