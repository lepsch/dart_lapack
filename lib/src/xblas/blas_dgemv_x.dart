// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/matrix.dart';

void Function(
  int trans,
  int m,
  int n,
  double alpha,
  Matrix<double> A,
  int lda,
  Array<double> x,
  int incx,
  double beta,
  Array<double> y,
  int incy,
  int prec,
) blas_dgemv_x = (
  trans,
  m,
  n,
  alpha,
  A,
  lda,
  x,
  incx,
  beta,
  y,
  incy,
  prec,
) {};
