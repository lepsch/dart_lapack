// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlaswp.dart';

void zgesc2(
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> RHS_,
  final Array<int> IPIV_,
  final Array<int> JPIV_,
  final Box<double> SCALE,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final RHS = RHS_.having();
  final IPIV = IPIV_.having();
  final JPIV = JPIV_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  int I, J;
  double EPS, SMLNUM;
  Complex TEMP;

  // Set constant to control overflow

  EPS = dlamch('P');
  SMLNUM = dlamch('S') / EPS;

  // Apply permutations IPIV to RHS

  zlaswp(1, RHS.asMatrix(LDA), LDA, 1, N - 1, IPIV, 1);

  // Solve for L part

  for (I = 1; I <= N - 1; I++) {
    for (J = I + 1; J <= N; J++) {
      RHS[J] -= A[J][I] * RHS[I];
    }
  }

  // Solve for U part

  SCALE.value = ONE;

  // Check for scaling

  I = izamax(N, RHS, 1);
  if (TWO * SMLNUM * RHS[I].abs() > A[N][N].abs()) {
    TEMP = Complex(ONE / TWO, ZERO) / RHS[I].abs().toComplex();
    zscal(N, TEMP, RHS(1), 1);
    SCALE.value *= TEMP.real;
  }
  for (I = N; I >= 1; I--) {
    TEMP = Complex(ONE, ZERO) / A[I][I];
    RHS[I] *= TEMP;
    for (J = I + 1; J <= N; J++) {
      RHS[I] -= RHS[J] * (A[I][J] * TEMP);
    }
  }

  // Apply permutations JPIV to the solution (RHS)

  zlaswp(1, RHS.asMatrix(LDA), LDA, 1, N - 1, JPIV, -1);
}
