// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/dtrmm.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dlarzb(
  final String SIDE,
  final String TRANS,
  final String DIRECT,
  final String STOREV,
  final int M,
  final int N,
  final int K,
  final int L,
  final Matrix<double> V_,
  final int LDV,
  final Matrix<double> T_,
  final int LDT,
  final Matrix<double> C_,
  final int LDC,
  final Matrix<double> WORK_,
  final int LDWORK,
) {
  final V = V_.having(ld: LDV);
  final T = T_.having(ld: LDT);
  final C = C_.having(ld: LDC);
  final WORK = WORK_.having(ld: LDWORK);

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0;
  String TRANST;
  int I, INFO, J;

  // Quick return if possible

  if (M <= 0 || N <= 0) return;

  // Check for currently supported options

  INFO = 0;
  if (!lsame(DIRECT, 'B')) {
    INFO = -3;
  } else if (!lsame(STOREV, 'R')) {
    INFO = -4;
  }
  if (INFO != 0) {
    xerbla('DLARZB', -INFO);
    return;
  }

  if (lsame(TRANS, 'N')) {
    TRANST = 'T';
  } else {
    TRANST = 'N';
  }

  if (lsame(SIDE, 'L')) {
    // Form  H * C  or  H**T * C

    // W( 1:n, 1:k ) = C( 1:k, 1:n )**T

    for (J = 1; J <= K; J++) {
      dcopy(N, C(J, 1).asArray(), LDC, WORK(1, J).asArray(), 1);
    }

    // W( 1:n, 1:k ) = W( 1:n, 1:k ) + ...
    //                 C( m-l+1:m, 1:n )**T * V( 1:k, 1:l )**T

    if (L > 0) {
      dgemm('Transpose', 'Transpose', N, K, L, ONE, C(M - L + 1, 1), LDC, V,
          LDV, ONE, WORK, LDWORK);
    }

    // W( 1:n, 1:k ) = W( 1:n, 1:k ) * T**T  or  W( 1:m, 1:k ) * T

    dtrmm(
        'Right', 'Lower', TRANST, 'Non-unit', N, K, ONE, T, LDT, WORK, LDWORK);

    // C( 1:k, 1:n ) = C( 1:k, 1:n ) - W( 1:n, 1:k )**T

    for (J = 1; J <= N; J++) {
      for (I = 1; I <= K; I++) {
        C[I][J] -= WORK[J][I];
      }
    }

    // C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
    //                     V( 1:k, 1:l )**T * W( 1:n, 1:k )**T

    if (L > 0) {
      dgemm('Transpose', 'Transpose', L, N, K, -ONE, V, LDV, WORK, LDWORK, ONE,
          C(M - L + 1, 1), LDC);
    }
  } else if (lsame(SIDE, 'R')) {
    // Form  C * H  or  C * H**T

    // W( 1:m, 1:k ) = C( 1:m, 1:k )

    for (J = 1; J <= K; J++) {
      dcopy(M, C(1, J).asArray(), 1, WORK(1, J).asArray(), 1);
    }

    // W( 1:m, 1:k ) = W( 1:m, 1:k ) + ...
    //                 C( 1:m, n-l+1:n ) * V( 1:k, 1:l )**T

    if (L > 0) {
      dgemm('No transpose', 'Transpose', M, K, L, ONE, C(1, N - L + 1), LDC, V,
          LDV, ONE, WORK, LDWORK);
    }

    // W( 1:m, 1:k ) = W( 1:m, 1:k ) * T  or  W( 1:m, 1:k ) * T**T

    dtrmm('Right', 'Lower', TRANS, 'Non-unit', M, K, ONE, T, LDT, WORK, LDWORK);

    // C( 1:m, 1:k ) = C( 1:m, 1:k ) - W( 1:m, 1:k )

    for (J = 1; J <= K; J++) {
      for (I = 1; I <= M; I++) {
        C[I][J] -= WORK[I][J];
      }
    }

    // C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
    //                     W( 1:m, 1:k ) * V( 1:k, 1:l )

    if (L > 0) {
      dgemm('No transpose', 'No transpose', M, L, K, -ONE, WORK, LDWORK, V, LDV,
          ONE, C(1, N - L + 1), LDC);
    }
  }
}
