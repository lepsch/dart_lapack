// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';

void zget54(
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> S_,
  final int LDS,
  final Matrix<Complex> T_,
  final int LDT,
  final Matrix<Complex> U_,
  final int LDU,
  final Matrix<Complex> V_,
  final int LDV,
  final Array<Complex> WORK_,
  final Box<double> RESULT,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final S = S_.having(ld: LDS);
  final T = T_.having(ld: LDT);
  final U = U_.having(ld: LDU);
  final V = V_.having(ld: LDV);
  final WORK = WORK_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  double ABNORM, ULP, UNFL, WNORM;
  final DUM = Array<double>(1);

  RESULT.value = ZERO;
  if (N <= 0) return;

  // Constants

  UNFL = dlamch('Safe minimum');
  ULP = dlamch('Epsilon') * dlamch('Base');

  // compute the norm of (A,B)

  zlacpy('Full', N, N, A, LDA, WORK.asMatrix(), N);
  zlacpy('Full', N, N, B, LDB, WORK(N * N + 1).asMatrix(), N);
  ABNORM = max(zlange('1', N, 2 * N, WORK.asMatrix(), N, DUM), UNFL);

  // Compute W1 = A - U*S*V', and put in the array WORK(1:N*N)

  zlacpy(' ', N, N, A, LDA, WORK.asMatrix(), N);
  zgemm('N', 'N', N, N, N, Complex.one, U, LDU, S, LDS, Complex.zero,
      WORK(N * N + 1).asMatrix(), N);

  zgemm('N', 'C', N, N, N, -Complex.one, WORK(N * N + 1).asMatrix(), N, V, LDV,
      Complex.one, WORK.asMatrix(), N);

  // Compute W2 = B - U*T*V', and put in the workarray W(N*N+1:2*N*N)

  zlacpy(' ', N, N, B, LDB, WORK(N * N + 1).asMatrix(), N);
  zgemm('N', 'N', N, N, N, Complex.one, U, LDU, T, LDT, Complex.zero,
      WORK(2 * N * N + 1).asMatrix(), N);

  zgemm('N', 'C', N, N, N, -Complex.one, WORK(2 * N * N + 1).asMatrix(), N, V,
      LDV, Complex.one, WORK(N * N + 1).asMatrix(), N);

  // Compute norm(W)/ ( ulp*norm((A,B)) )

  WNORM = zlange('1', N, 2 * N, WORK.asMatrix(), N, DUM);

  if (ABNORM > WNORM) {
    RESULT.value = (WNORM / ABNORM) / (2 * N * ULP);
  } else {
    if (ABNORM < ONE) {
      RESULT.value = (min(WNORM, 2 * N * ABNORM) / ABNORM) / (2 * N * ULP);
    } else {
      RESULT.value = min(WNORM / ABNORM, (2 * N)) / (2 * N * ULP);
    }
  }
}
