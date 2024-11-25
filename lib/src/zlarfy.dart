// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/zaxpy.dart';
import 'package:dart_lapack/src/blas/zdotc.dart';
import 'package:dart_lapack/src/blas/zhemv.dart';
import 'package:dart_lapack/src/blas/zher2.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

void zlarfy(
  final String UPLO,
  final int N,
  final Array<Complex> V_,
  final int INCV,
  final Complex TAU,
  final Matrix<Complex> C_,
  final int LDC,
  final Array<Complex> WORK_,
) {
  final C = C_.having(ld: LDC);
  final V = V_.having();
  final WORK = WORK_.having();
  const HALF = Complex(0.5, 0.0);
  Complex ALPHA;

  if (TAU == Complex.zero) return;

  // Form  w:= C * v

  zhemv(UPLO, N, Complex.one, C, LDC, V, INCV, Complex.zero, WORK, 1);

  ALPHA = -HALF * TAU * zdotc(N, WORK, 1, V, INCV);
  zaxpy(N, ALPHA, V, INCV, WORK, 1);

  // C := C - v * w' - w * v'

  zher2(UPLO, N, -TAU, V, INCV, WORK, 1, C, LDC);
}
