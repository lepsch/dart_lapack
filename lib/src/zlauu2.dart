// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zdotc.dart';
import 'package:dart_lapack/src/blas/zdscal.dart';
import 'package:dart_lapack/src/blas/zgemv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlacgv.dart';

void zlauu2(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  bool UPPER;
  int I;
  double AII;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('ZLAUU2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (UPPER) {
    // Compute the product U * U**H.

    for (I = 1; I <= N; I++) {
      AII = A[I][I].real;
      if (I < N) {
        A[I][I] = (AII * AII +
                zdotc(N - I, A(I, I + 1).asArray(), LDA, A(I, I + 1).asArray(),
                        LDA)
                    .real)
            .toComplex();
        zlacgv(N - I, A(I, I + 1).asArray(), LDA);
        zgemv('No transpose', I - 1, N - I, Complex.one, A(1, I + 1), LDA,
            A(I, I + 1).asArray(), LDA, AII.toComplex(), A(1, I).asArray(), 1);
        zlacgv(N - I, A(I, I + 1).asArray(), LDA);
      } else {
        zdscal(I, AII, A(1, I).asArray(), 1);
      }
    }
  } else {
    // Compute the product L**H * L.

    for (I = 1; I <= N; I++) {
      AII = A[I][I].real;
      if (I < N) {
        A[I][I] = (AII * AII +
                zdotc(N - I, A(I + 1, I).asArray(), 1, A(I + 1, I).asArray(), 1)
                    .real)
            .toComplex();
        zlacgv(I - 1, A(I, 1).asArray(), LDA);
        zgemv(
            'Conjugate transpose',
            N - I,
            I - 1,
            Complex.one,
            A(I + 1, 1),
            LDA,
            A(I + 1, I).asArray(),
            1,
            AII.toComplex(),
            A(I, 1).asArray(),
            LDA);
        zlacgv(I - 1, A(I, 1).asArray(), LDA);
      } else {
        zdscal(I, AII, A(I, 1).asArray(), LDA);
      }
    }
  }
}
