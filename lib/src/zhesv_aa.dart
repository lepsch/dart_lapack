// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhetrf_aa.dart';
import 'package:lapack/src/zhetrs_aa.dart';

void zhesv_aa(
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
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  final WORK = WORK_.having();

  bool LQUERY;
  int LWKMIN, LWKOPT = 0, LWKOPT_HETRF, LWKOPT_HETRS;

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
    zhetrf_aa(UPLO, N, A, LDA, IPIV, WORK, -1, INFO);
    LWKOPT_HETRF = WORK[1].toInt();
    zhetrs_aa(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, -1, INFO);
    LWKOPT_HETRS = WORK[1].toInt();
    LWKOPT = max(LWKMIN, max(LWKOPT_HETRF, LWKOPT_HETRS));
    WORK[1] = LWKOPT.toComplex();
  }

  if (INFO.value != 0) {
    xerbla('ZHESV_AA', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Compute the factorization A = U**H*T*U or A = L*T*L**H.

  zhetrf_aa(UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO);
  if (INFO.value == 0) {
    // Solve the system A*X = B, overwriting B with X.

    zhetrs_aa(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO);
  }

  WORK[1] = LWKOPT.toComplex();
}
