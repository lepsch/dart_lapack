// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacn2.dart';
import 'package:dart_lapack/src/dsytrs_3.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dsycon_3(
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> E_,
  final Array<int> IPIV_,
  final double ANORM,
  final Box<double> RCOND,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final E = E_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  bool UPPER;
  int I;
  final KASE = Box(0);
  final AINVNM = Box(0.0);
  final ISAVE = Array<int>(3);

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  } else if (ANORM < ZERO) {
    INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('DSYCON_3', -INFO.value);
    return;
  }

  // Quick return if possible

  RCOND.value = ZERO;
  if (N == 0) {
    RCOND.value = ONE;
    return;
  } else if (ANORM <= ZERO) {
    return;
  }

  // Check that the diagonal matrix D is nonsingular.

  if (UPPER) {
    // Upper triangular storage: examine D from bottom to top

    for (I = N; I >= 1; I--) {
      if (IPIV[I] > 0 && A[I][I] == ZERO) return;
    }
  } else {
    // Lower triangular storage: examine D from top to bottom.

    for (I = 1; I <= N; I++) {
      if (IPIV[I] > 0 && A[I][I] == ZERO) return;
    }
  }

  // Estimate the 1-norm of the inverse.

  KASE.value = 0;
  while (true) {
    dlacn2(N, WORK(N + 1), WORK, IWORK, AINVNM, KASE, ISAVE);
    if (KASE.value == 0) break;

    // Multiply by inv(L*D*L**T) or inv(U*D*U**T).

    dsytrs_3(UPLO, N, 1, A, LDA, E, IPIV, WORK.asMatrix(N), N, INFO);
  }

  // Compute the estimate of the reciprocal condition number.

  if (AINVNM.value != ZERO) RCOND.value = (ONE / AINVNM.value) / ANORM;
}
