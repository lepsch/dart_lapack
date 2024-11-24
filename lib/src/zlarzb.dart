// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/ztrmm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacgv.dart';

void zlarzb(
  final String SIDE,
  final String TRANS,
  final String DIRECT,
  final String STOREV,
  final int M,
  final int N,
  final int K,
  final int L,
  final Matrix<Complex> V_,
  final int LDV,
  final Matrix<Complex> T_,
  final int LDT,
  final Matrix<Complex> C_,
  final int LDC,
  final Matrix<Complex> WORK_,
  final int LDWORK,
) {
  final V = V_.having(ld: LDV);
  final T = T_.having(ld: LDT);
  final C = C_.having(ld: LDC);
  final WORK = WORK_.having(ld: LDWORK);

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  String TRANST;
  int I, J;
  final INFO = Box(0);

  // ..
  // .. External Functions ..
  //- bool               lsame;
  // EXTERNAL lsame
  // ..
  // .. External Subroutines ..
  // EXTERNAL XERBLA, ZCOPY, ZGEMM, ZLACGV, ZTRMM

  // Quick return if possible

  if (M <= 0 || N <= 0) return;

  // Check for currently supported options

  INFO.value = 0;
  if (!lsame(DIRECT, 'B')) {
    INFO.value = -3;
  } else if (!lsame(STOREV, 'R')) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('ZLARZB', -INFO.value);
    return;
  }

  if (lsame(TRANS, 'N')) {
    TRANST = 'C';
  } else {
    TRANST = 'N';
  }

  if (lsame(SIDE, 'L')) {
    // Form  H * C  or  H**H * C

    // W( 1:n, 1:k ) = C( 1:k, 1:n )**H

    for (J = 1; J <= K; J++) {
      zcopy(N, C(J, 1).asArray(), LDC, WORK(1, J).asArray(), 1);
    }

    // W( 1:n, 1:k ) = W( 1:n, 1:k ) + ...
    //                 C( m-l+1:m, 1:n )**H * V( 1:k, 1:l )**T

    if (L > 0) {
      zgemm('Transpose', 'Conjugate transpose', N, K, L, Complex.one,
          C(M - L + 1, 1), LDC, V, LDV, Complex.one, WORK, LDWORK);
    }

    // W( 1:n, 1:k ) = W( 1:n, 1:k ) * T**T  or  W( 1:m, 1:k ) * T

    ztrmm('Right', 'Lower', TRANST, 'Non-unit', N, K, Complex.one, T, LDT, WORK,
        LDWORK);

    // C( 1:k, 1:n ) = C( 1:k, 1:n ) - W( 1:n, 1:k )**H

    for (J = 1; J <= N; J++) {
      for (I = 1; I <= K; I++) {
        C[I][J] -= WORK[J][I];
      }
    }

    // C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
    //                     V( 1:k, 1:l )**H * W( 1:n, 1:k )**H

    if (L > 0) {
      zgemm('Transpose', 'Transpose', L, N, K, -Complex.one, V, LDV, WORK,
          LDWORK, Complex.one, C(M - L + 1, 1), LDC);
    }
  } else if (lsame(SIDE, 'R')) {
    // Form  C * H  or  C * H**H

    // W( 1:m, 1:k ) = C( 1:m, 1:k )

    for (J = 1; J <= K; J++) {
      zcopy(M, C(1, J).asArray(), 1, WORK(1, J).asArray(), 1);
    }

    // W( 1:m, 1:k ) = W( 1:m, 1:k ) + ...
    //                 C( 1:m, n-l+1:n ) * V( 1:k, 1:l )**H

    if (L > 0) {
      zgemm('No transpose', 'Transpose', M, K, L, Complex.one, C(1, N - L + 1),
          LDC, V, LDV, Complex.one, WORK, LDWORK);
    }

    // W( 1:m, 1:k ) = W( 1:m, 1:k ) * conjg( T )  or
    //                 W( 1:m, 1:k ) * T**H

    for (J = 1; J <= K; J++) {
      zlacgv(K - J + 1, T(J, J).asArray(), 1);
    }
    ztrmm('Right', 'Lower', TRANS, 'Non-unit', M, K, Complex.one, T, LDT, WORK,
        LDWORK);
    for (J = 1; J <= K; J++) {
      zlacgv(K - J + 1, T(J, J).asArray(), 1);
    }

    // C( 1:m, 1:k ) = C( 1:m, 1:k ) - W( 1:m, 1:k )

    for (J = 1; J <= K; J++) {
      for (I = 1; I <= M; I++) {
        C[I][J] -= WORK[I][J];
      }
    }

    // C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
    //                     W( 1:m, 1:k ) * conjg( V( 1:k, 1:l ) )

    for (J = 1; J <= L; J++) {
      zlacgv(K, V(1, J).asArray(), 1);
    }
    if (L > 0) {
      zgemm('No transpose', 'No transpose', M, L, K, -Complex.one, WORK, LDWORK,
          V, LDV, Complex.one, C(1, N - L + 1), LDC);
    }
    for (J = 1; J <= L; J++) {
      zlacgv(K, V(1, J).asArray(), 1);
    }
  }
}
