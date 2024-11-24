// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/matrix.dart';

void Function(
  int trans,
  int m,
  int n,
  double alpha,
  Matrix<double> A,
  int lda,
  Array<double> head_x,
  Array<double> tail_x,
  int incx,
  double beta,
  Array<double> y,
  int incy,
  int prec,
) blas_dgemv2_x = (
  trans,
  m,
  n,
  alpha,
  A,
  lda,
  head_x,
  tail_x,
  incx,
  beta,
  y,
  incy,
  prec,
) {};
