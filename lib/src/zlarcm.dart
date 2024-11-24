// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zlarcm(
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> C_,
  final int LDC,
  final Array<double> RWORK_,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final C = C_.having(ld: LDC);
  final RWORK = RWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  int I, J, L;

  // Quick return if possible.

  if ((M == 0) || (N == 0)) return;

  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      RWORK[(J - 1) * M + I] = B[I][J].real;
    }
  }

  L = M * N + 1;
  dgemm('N', 'N', M, N, M, ONE, A, LDA, RWORK.asMatrix(M), M, ZERO,
      RWORK(L).asMatrix(M), M);
  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      C[I][J] = RWORK[L + (J - 1) * M + I - 1].toComplex();
    }
  }

  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      RWORK[(J - 1) * M + I] = B[I][J].imaginary;
    }
  }
  dgemm('N', 'N', M, N, M, ONE, A, LDA, RWORK.asMatrix(M), M, ZERO,
      RWORK(L).asMatrix(M), M);
  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      C[I][J] = Complex(C[I][J].real, RWORK[L + (J - 1) * M + I - 1]);
    }
  }
}
