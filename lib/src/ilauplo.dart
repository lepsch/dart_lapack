// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/install/lsame.dart';

int ilauplo(final String UPLO) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const BLAS_UPPER = 121, BLAS_LOWER = 122;

  if (lsame(UPLO, 'U')) {
    return BLAS_UPPER;
  } else if (lsame(UPLO, 'L')) {
    return BLAS_LOWER;
  } else {
    return -1;
  }
}
