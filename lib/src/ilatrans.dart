// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/install/lsame.dart';

int ilatrans(final String TRANS) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const BLAS_NO_TRANS = 111, BLAS_TRANS = 112, BLAS_CONJ_TRANS = 113;

  if (lsame(TRANS, 'N')) {
    return BLAS_NO_TRANS;
  } else if (lsame(TRANS, 'T')) {
    return BLAS_TRANS;
  } else if (lsame(TRANS, 'C')) {
    return BLAS_CONJ_TRANS;
  } else {
    return -1;
  }
}
