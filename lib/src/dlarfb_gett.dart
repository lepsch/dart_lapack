// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/dtrmm.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';

void dlarfb_gett(
  final String IDENT,
  final int M,
  final int N,
  final int K,
  final Matrix<double> T_,
  final int LDT,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> WORK_,
  final int LDWORK,
) {
  final T = T_.having(ld: LDT);
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final WORK = WORK_.having(ld: LDWORK);
  const ONE = 1.0, ZERO = 0.0;
  bool LNOTIDENT;
  int I, J;

  // Quick return if possible
  if (M < 0 || N <= 0 || K == 0 || K > N) return;

  LNOTIDENT = !lsame(IDENT, 'I');

  // First Step. Computation of the Column Block 2:
  //
  // ( A2 ) := H * ( A2 )
  // ( B2 )        ( B2 )

  if (N > K) {
    // col2_(1) Compute W2: = A2. Therefore, copy A2 = A(1:K, K+1:N)
    // into W2=WORK(1:K, 1:N-K) column-by-column.
    for (J = 1; J <= N - K; J++) {
      dcopy(K, A(1, K + J).asArray(), 1, WORK(1, J).asArray(), 1);
    }

    if (LNOTIDENT) {
      // col2_(2) Compute W2: = (V1**T) * W2 = (A1**T) * W2,
      // V1 is not an identity matrix, but unit lower-triangular
      // V1 stored in A1 (diagonal ones are not stored).
      dtrmm('L', 'L', 'T', 'U', K, N - K, ONE, A, LDA, WORK, LDWORK);
    }

    // col2_(3) Compute W2: = W2 + (V2**T) * B2 = W2 + (B1**T) * B2
    // V2 stored in B1.
    if (M > 0) {
      dgemm('T', 'N', K, N - K, M, ONE, B, LDB, B(1, K + 1), LDB, ONE, WORK,
          LDWORK);
    }

    // col2_(4) Compute W2: = T * W2,
    // T is upper-triangular.
    dtrmm('L', 'U', 'N', 'N', K, N - K, ONE, T, LDT, WORK, LDWORK);

    // col2_(5) Compute B2: = B2 - V2 * W2 = B2 - B1 * W2,
    // V2 stored in B1.
    if (M > 0) {
      dgemm('N', 'N', M, N - K, K, -ONE, B, LDB, WORK, LDWORK, ONE, B(1, K + 1),
          LDB);
    }

    if (LNOTIDENT) {
      // col2_(6) Compute W2: = V1 * W2 = A1 * W2,
      // V1 is not an identity matrix, but unit lower-triangular,
      // V1 stored in A1 (diagonal ones are not stored).
      dtrmm('L', 'L', 'N', 'U', K, N - K, ONE, A, LDA, WORK, LDWORK);
    }

    // col2_(7) Compute A2: = A2 - W2 =
    //                      = A(1:K, K+1:N-K) - WORK(1:K, 1:N-K),
    // column-by-column.
    for (J = 1; J <= N - K; J++) {
      for (I = 1; I <= K; I++) {
        A[I][K + J] -= WORK[I][J];
      }
    }
  }

  // Second Step. Computation of the Column Block 1:
  //
  // ( A1 ) := H * ( A1 )
  // ( B1 )        (  0 )

  // col1_(1) Compute W1: = A1. Copy the upper-triangular
  // A1 = A(1:K, 1:K) into the upper-triangular
  // W1 = WORK(1:K, 1:K) column-by-column.
  for (J = 1; J <= K; J++) {
    dcopy(J, A(1, J).asArray(), 1, WORK(1, J).asArray(), 1);
  }

  // Set the subdiagonal elements of W1 to zero column-by-column.
  for (J = 1; J <= K - 1; J++) {
    for (I = J + 1; I <= K; I++) {
      WORK[I][J] = ZERO;
    }
  }

  if (LNOTIDENT) {
    // col1_(2) Compute W1: = (V1**T) * W1 = (A1**T) * W1,
    // V1 is not an identity matrix, but unit lower-triangular
    // V1 stored in A1 (diagonal ones are not stored),
    // W1 is upper-triangular with zeroes below the diagonal.
    dtrmm('L', 'L', 'T', 'U', K, K, ONE, A, LDA, WORK, LDWORK);
  }

  // col1_(3) Compute W1: = T * W1,
  // T is upper-triangular,
  // W1 is upper-triangular with zeroes below the diagonal.
  dtrmm('L', 'U', 'N', 'N', K, K, ONE, T, LDT, WORK, LDWORK);

  // col1_(4) Compute B1: = - V2 * W1 = - B1 * W1,
  // V2 = B1, W1 is upper-triangular with zeroes below the diagonal.
  if (M > 0) {
    dtrmm('R', 'U', 'N', 'N', M, K, -ONE, WORK, LDWORK, B, LDB);
  }

  if (LNOTIDENT) {
    // col1_(5) Compute W1: = V1 * W1 = A1 * W1,
    // V1 is not an identity matrix, but unit lower-triangular
    // V1 stored in A1 (diagonal ones are not stored),
    // W1 is upper-triangular on input with zeroes below the diagonal,
    // and square on output.
    dtrmm('L', 'L', 'N', 'U', K, K, ONE, A, LDA, WORK, LDWORK);

    // col1_(6) Compute A1: = A1 - W1 = A(1:K, 1:K) - WORK(1:K, 1:K)
    // column-by-column. A1 is upper-triangular on input.
    // If IDENT, A1 is square on output, and W1 is square,
    // if NOT IDENT, A1 is upper-triangular on output,
    // W1 is upper-triangular.

    // col1_(6)_a Compute elements of A1 below the diagonal.
    for (J = 1; J <= K - 1; J++) {
      for (I = J + 1; I <= K; I++) {
        A[I][J] = -WORK[I][J];
      }
    }
  }

  // col1_(6)_b Compute elements of A1 on and above the diagonal.
  for (J = 1; J <= K; J++) {
    for (I = 1; I <= J; I++) {
      A[I][J] -= WORK[I][J];
    }
  }
}
