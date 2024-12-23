// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dznrm2.dart';
import 'package:dart_lapack/src/blas/zgemv.dart';
import 'package:dart_lapack/src/blas/zgerc.dart';
import 'package:dart_lapack/src/blas/zscal.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlarnv.dart';

void zlarge(
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> ISEED_,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final ISEED = ISEED_.having(length: 4);
  final WORK = WORK_.having();
  const ZERO = 0.0;
  int I;
  double WN;
  Complex TAU, WA, WB;

  // Test the input arguments

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
  } else if (LDA < max(1, N)) {
    INFO.value = -3;
  }
  if (INFO.value < 0) {
    xerbla('ZLARGE', -INFO.value);
    return;
  }

  // pre- and post-multiply A by random unitary matrix

  for (I = N; I >= 1; I--) {
    // generate random reflection

    zlarnv(3, ISEED, N - I + 1, WORK);
    WN = dznrm2(N - I + 1, WORK, 1);
    WA = (WN / WORK[1].abs()).toComplex() * WORK[1];
    if (WN == ZERO) {
      TAU = Complex.zero;
    } else {
      WB = WORK[1] + WA;
      zscal(N - I, Complex.one / WB, WORK(2), 1);
      WORK[1] = Complex.one;
      TAU = (WB / WA).real.toComplex();
    }

    // multiply A(i:n,1:n) by random reflection from the left

    zgemv('Conjugate transpose', N - I + 1, N, Complex.one, A(I, 1), LDA, WORK,
        1, Complex.zero, WORK(N + 1), 1);
    zgerc(N - I + 1, N, -TAU, WORK, 1, WORK(N + 1), 1, A(I, 1), LDA);

    // multiply A(1:n,i:n) by random reflection from the right

    zgemv('No transpose', N, N - I + 1, Complex.one, A(1, I), LDA, WORK, 1,
        Complex.zero, WORK(N + 1), 1);
    zgerc(N, N - I + 1, -TAU, WORK(N + 1), 1, WORK, 1, A(1, I), LDA);
  }
}
