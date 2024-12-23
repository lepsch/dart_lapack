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

void ztpqrt2(
  final int M,
  final int N,
  final int L,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> T_,
  final int LDT,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final T = T_.having(ld: LDT);
  int I, J, P, MP, NP;
  Complex ALPHA;

  // Test the input arguments

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (L < 0 || L > min(M, N)) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDB < max(1, M)) {
    INFO.value = -7;
  } else if (LDT < max(1, N)) {
    INFO.value = -9;
  }
  if (INFO.value != 0) {
    xerbla('ZTPQRT2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || M == 0) return;

  for (I = 1; I <= N; I++) {
    // Generate elementary reflector H(I) to annihilate B(:,I)

    P = M - L + min(L, I);
    zlarfg(P + 1, A(I, I), B(1, I).asArray(), 1, T(I, 1));
    if (I < N) {
      // W(1:N-I) := C(I:M,I+1:N)**H * C(I:M,I) [use W = T(:,N)]

      for (J = 1; J <= N - I; J++) {
        T[J][N] = A[I][I + J].conjugate();
      }
      zgemv('C', P, N - I, Complex.one, B(1, I + 1), LDB, B(1, I).asArray(), 1,
          Complex.one, T(1, N).asArray(), 1);

      // C(I:M,I+1:N) = C(I:m,I+1:N) + alpha*C(I:M,I)*W(1:N-1)**H

      ALPHA = -T[I][1].conjugate();
      for (J = 1; J <= N - I; J++) {
        A[I][I + J] += ALPHA * T[J][N].conjugate();
      }
      zgerc(P, N - I, ALPHA, B(1, I).asArray(), 1, T(1, N).asArray(), 1,
          B(1, I + 1), LDB);
    }
  }

  for (I = 2; I <= N; I++) {
    // T(1:I-1,I) := C(I:M,1:I-1)**H * (alpha * C(I:M,I))

    ALPHA = -T[I][1];

    for (J = 1; J <= I - 1; J++) {
      T[J][I] = Complex.zero;
    }
    P = min(I - 1, L);
    MP = min(M - L + 1, M);
    NP = min(P + 1, N);

    // Triangular part of B2

    for (J = 1; J <= P; J++) {
      T[J][I] = ALPHA * B[M - L + J][I];
    }
    ztrmv('U', 'C', 'N', P, B(MP, 1), LDB, T(1, I).asArray(), 1);

    // Rectangular part of B2

    zgemv('C', L, I - 1 - P, ALPHA, B(MP, NP), LDB, B(MP, I).asArray(), 1,
        Complex.zero, T(NP, I).asArray(), 1);

    // B1

    zgemv('C', M - L, I - 1, ALPHA, B, LDB, B(1, I).asArray(), 1, Complex.one,
        T(1, I).asArray(), 1);

    // T(1:I-1,I) := T(1:I-1,1:I-1) * T(1:I-1,I)

    ztrmv('U', 'N', 'N', I - 1, T, LDT, T(1, I).asArray(), 1);

    // T(I,I) = tau(I)

    T[I][I] = T[I][1];
    T[I][1] = Complex.zero;
  }
}
