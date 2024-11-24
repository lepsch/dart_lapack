// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zgemv.dart';
import 'package:dart_lapack/src/blas/zgerc.dart';
import 'package:dart_lapack/src/blas/ztrmv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlarfg.dart';

void zgeqrt2(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> T_,
  final int LDT,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final T = T_.having(ld: LDT);
  int I, K;
  Complex AII, ALPHA;

  // Test the input arguments

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -2;
  } else if (M < N) {
    INFO.value = -1;
  } else if (LDA < max(1, M)) {
    INFO.value = -4;
  } else if (LDT < max(1, N)) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('ZGEQRT2', -INFO.value);
    return;
  }

  K = min(M, N);

  for (I = 1; I <= K; I++) {
    // Generate elem. refl. H(i) to annihilate A(i+1:m,i), tau(I) -> T(I,1)

    zlarfg(M - I + 1, A(I, I), A(min(I + 1, M), I).asArray(), 1, T(I, 1));
    if (I < N) {
      // Apply H(i) to A(I:M,I+1:N) from the left

      AII = A[I][I];
      A[I][I] = Complex.one;

      // W(1:N-I) := A(I:M,I+1:N)^H * A(I:M,I) [W = T(:,N)]

      zgemv('C', M - I + 1, N - I, Complex.one, A(I, I + 1), LDA,
          A(I, I).asArray(), 1, Complex.zero, T(1, N).asArray(), 1);

      // A(I:M,I+1:N) = A(I:m,I+1:N) + alpha*A(I:M,I)*W(1:N-1)^H

      ALPHA = -T[I][1].conjugate();
      zgerc(M - I + 1, N - I, ALPHA, A(I, I).asArray(), 1, T(1, N).asArray(), 1,
          A(I, I + 1), LDA);
      A[I][I] = AII;
    }
  }

  for (I = 2; I <= N; I++) {
    AII = A[I][I];
    A[I][I] = Complex.one;

    // T(1:I-1,I) := alpha * A(I:M,1:I-1)**H * A(I:M,I)

    ALPHA = -T[I][1];
    zgemv('C', M - I + 1, I - 1, ALPHA, A(I, 1), LDA, A(I, I).asArray(), 1,
        Complex.zero, T(1, I).asArray(), 1);
    A[I][I] = AII;

    // T(1:I-1,I) := T(1:I-1,1:I-1) * T(1:I-1,I)

    ztrmv('U', 'N', 'N', I - 1, T, LDT, T(1, I).asArray(), 1);

    // T(I,I) = tau(I)

    T[I][I] = T[I][1];
    T[I][1] = Complex.zero;
  }
}
