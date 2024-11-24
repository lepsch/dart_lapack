// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zgeru.dart';
import 'package:dart_lapack/src/blas/zswap.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void zgetc2(
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Array<int> JPIV_,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final JPIV = JPIV_.having();

  const ZERO = 0.0, ONE = 1.0;
  int I, IP, IPV = 0, J, JP, JPV = 0;
  double EPS, SMIN = 0, SMLNUM, XMAX;

  INFO.value = 0;

  // Quick return if possible

  if (N == 0) return;

  // Set constants to control overflow

  EPS = dlamch('P');
  SMLNUM = dlamch('S') / EPS;

  // Handle the case N=1 by itself

  if (N == 1) {
    IPIV[1] = 1;
    JPIV[1] = 1;
    if (A[1][1].abs() < SMLNUM) {
      INFO.value = 1;
      A[1][1] = Complex(SMLNUM, ZERO);
    }
    return;
  }

  // Factorize A using complete pivoting.
  // Set pivots less than SMIN to SMIN

  for (I = 1; I <= N - 1; I++) {
    // Find max element in matrix A

    XMAX = ZERO;
    for (IP = I; IP <= N; IP++) {
      for (JP = I; JP <= N; JP++) {
        if (A[IP][JP].abs() >= XMAX) {
          XMAX = A[IP][JP].abs();
          IPV = IP;
          JPV = JP;
        }
      }
    }
    if (I == 1) SMIN = max(EPS * XMAX, SMLNUM);

    // Swap rows

    if (IPV != I) zswap(N, A(IPV, 1).asArray(), LDA, A(I, 1).asArray(), LDA);
    IPIV[I] = IPV;

    // Swap columns

    if (JPV != I) zswap(N, A(1, JPV).asArray(), 1, A(1, I).asArray(), 1);
    JPIV[I] = JPV;

    // Check for singularity

    if (A[I][I].abs() < SMIN) {
      INFO.value = I;
      A[I][I] = Complex(SMIN, ZERO);
    }
    for (J = I + 1; J <= N; J++) {
      A[J][I] /= A[I][I];
    }
    zgeru(N - I, N - I, -Complex(ONE), A(I + 1, I).asArray(), 1,
        A(I, I + 1).asArray(), LDA, A(I + 1, I + 1), LDA);
  }

  if (A[N][N].abs() < SMIN) {
    INFO.value = N;
    A[N][N] = Complex(SMIN, ZERO);
  }

  // Set last pivots to N

  IPIV[N] = N;
  JPIV[N] = N;
}
