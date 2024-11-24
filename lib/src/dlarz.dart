// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dger.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';

void dlarz(
  final String SIDE,
  final int M,
  final int N,
  final int L,
  final Array<double> V_,
  final int INCV,
  final double TAU,
  final Matrix<double> C_,
  final int LDC,
  final Array<double> WORK_,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final V = V_.having();
  final C = C_.having(ld: LDC);
  final WORK = WORK_.having();
  const ONE = 1.0, ZERO = 0.0;

  if (lsame(SIDE, 'L')) {
    // Form  H * C

    if (TAU != ZERO) {
      // w( 1:n ) = C( 1, 1:n )

      dcopy(N, C.asArray(), LDC, WORK, 1);

      // w( 1:n ) = w( 1:n ) + C( m-l+1:m, 1:n )**T * v( 1:l )

      dgemv(
          'Transpose', L, N, ONE, C(M - L + 1, 1), LDC, V, INCV, ONE, WORK, 1);

      // C( 1, 1:n ) = C( 1, 1:n ) - tau * w( 1:n )

      daxpy(N, -TAU, WORK, 1, C.asArray(), LDC);

      // C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
      //                     tau * v( 1:l ) * w( 1:n )**T

      dger(L, N, -TAU, V, INCV, WORK, 1, C(M - L + 1, 1), LDC);
    }
  } else {
    // Form  C * H

    if (TAU != ZERO) {
      // w( 1:m ) = C( 1:m, 1 )

      dcopy(M, C.asArray(), 1, WORK, 1);

      // w( 1:m ) = w( 1:m ) + C( 1:m, n-l+1:n, 1:n ) * v( 1:l )

      dgemv('No transpose', M, L, ONE, C(1, N - L + 1), LDC, V, INCV, ONE, WORK,
          1);

      // C( 1:m, 1 ) = C( 1:m, 1 ) - tau * w( 1:m )

      daxpy(M, -TAU, WORK, 1, C.asArray(), 1);

      // C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
      //                     tau * w( 1:m ) * v( 1:l )**T

      dger(M, L, -TAU, WORK, 1, V, INCV, C(1, N - L + 1), LDC);
    }
  }
}
