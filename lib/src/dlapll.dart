// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/daxpy.dart';
import 'package:dart_lapack/src/blas/ddot.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlarfg.dart';
import 'package:dart_lapack/src/dlas2.dart';
import 'package:dart_lapack/src/matrix.dart';

void dlapll(
  final int N,
  final Array<double> X_,
  final int INCX,
  final Array<double> Y_,
  final int INCY,
  final Box<double> SSMIN,
) {
  final X = X_.having();
  final Y = Y_.having();
  const ZERO = 0.0, ONE = 1.0;
  double A11, A12, A22, C;
  final SSMAX = Box(0.0), TAU = Box(0.0);

  // Quick return if possible

  if (N <= 1) {
    SSMIN.value = ZERO;
    return;
  }

  // Compute the QR factorization of the N-by-2 matrix ( X Y )

  dlarfg(N, X.box(1), X(1 + INCX), INCX, TAU);
  A11 = X[1];
  X[1] = ONE;

  C = -TAU.value * ddot(N, X, INCX, Y, INCY);
  daxpy(N, C, X, INCX, Y, INCY);

  dlarfg(N - 1, Y.box(1 + INCY), Y(1 + 2 * INCY), INCY, TAU);

  A12 = Y[1];
  A22 = Y[1 + INCY];

  // Compute the SVD of 2-by-2 Upper triangular matrix.

  dlas2(A11, A12, A22, SSMIN, SSMAX);
}
