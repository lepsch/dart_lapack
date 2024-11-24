// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/idamax.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void zptcon(
  final int N,
  final Array<double> D_,
  final Array<Complex> E_,
  final double ANORM,
  final Box<double> RCOND,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final E = E_.having();
  final RWORK = RWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  int I, IX;
  double AINVNM;

  // Test the input arguments.

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
  } else if (ANORM < ZERO) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('ZPTCON', -INFO.value);
    return;
  }

  // Quick return if possible

  RCOND.value = ZERO;
  if (N == 0) {
    RCOND.value = ONE;
    return;
  } else if (ANORM == ZERO) {
    return;
  }

  // Check that D(1:N) is positive.

  for (I = 1; I <= N; I++) {
    if (D[I] <= ZERO) return;
  }

  // Solve M(A) * x = e, where M(A) = (m(i,j)) is given by

  // m(i,j) =  abs(A(i,j)), i = j,
  // m(i,j) = -abs(A(i,j)), i != j,

  // and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**H.

  // Solve M(L) * x = e.

  RWORK[1] = ONE;
  for (I = 2; I <= N; I++) {
    RWORK[I] = ONE + RWORK[I - 1] * E[I - 1].abs();
  }

  // Solve D * M(L)**H * x = b.

  RWORK[N] /= D[N];
  for (I = N - 1; I >= 1; I--) {
    RWORK[I] = RWORK[I] / D[I] + RWORK[I + 1] * E[I].abs();
  }

  // Compute AINVNM = max(x[I]), 1<=i<=n.

  IX = idamax(N, RWORK, 1);
  AINVNM = RWORK[IX].abs();

  // Compute the reciprocal condition number.

  if (AINVNM != ZERO) RCOND.value = (ONE / AINVNM) / ANORM;
}
