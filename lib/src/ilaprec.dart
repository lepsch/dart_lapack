// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/install/lsame.dart';

int ilaprec(final String PREC) {
  const BLAS_PREC_SINGLE = 211,
      BLAS_PREC_DOUBLE = 212,
      BLAS_PREC_INDIGENOUS = 213,
      BLAS_PREC_EXTRA = 214;

  if (lsame(PREC, 'S')) {
    return BLAS_PREC_SINGLE;
  } else if (lsame(PREC, 'D')) {
    return BLAS_PREC_DOUBLE;
  } else if (lsame(PREC, 'I')) {
    return BLAS_PREC_INDIGENOUS;
  } else if (lsame(PREC, 'X') || lsame(PREC, 'E')) {
    return BLAS_PREC_EXTRA;
  } else {
    return -1;
  }
}
