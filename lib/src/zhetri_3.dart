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
import 'package:dart_lapack/src/zhetri_3x.dart';

void zhetri_3(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> E_,
  final Array<int> IPIV_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final E = E_.having();
  final IPIV = IPIV_.having();
  final WORK = WORK_.having();
  bool UPPER, LQUERY;
  int LWKOPT, NB;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  LQUERY = (LWORK == -1);

  // Determine the block size

  NB = max(1, ilaenv(1, 'ZHETRI_3', UPLO, N, -1, -1, -1));
  LWKOPT = (N + NB + 1) * (NB + 3);

  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  } else if (LWORK < LWKOPT && !LQUERY) {
    INFO.value = -8;
  }

  if (INFO.value != 0) {
    xerbla('ZHETRI_3', -INFO.value);
    return;
  } else if (LQUERY) {
    WORK[1] = LWKOPT.toComplex();
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  zhetri_3x(UPLO, N, A, LDA, E, IPIV, WORK.asMatrix(), NB, INFO);

  WORK[1] = LWKOPT.toComplex();
}
