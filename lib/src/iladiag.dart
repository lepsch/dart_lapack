// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/install/lsame.dart';

int iladiag(final String DIAG) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const BLAS_NON_UNIT_DIAG = 131, BLAS_UNIT_DIAG = 132;
  if (lsame(DIAG, 'N')) {
    return BLAS_NON_UNIT_DIAG;
  } else if (lsame(DIAG, 'U')) {
    return BLAS_UNIT_DIAG;
  } else {
    return -1;
  }
}
