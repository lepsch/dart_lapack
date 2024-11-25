// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zhptrs.dart';
import 'package:dart_lapack/src/zlacn2.dart';

void zhpcon(
  final String UPLO,
  final int N,
  final Array<Complex> AP,
  final Array<int> IPIV_,
  final double ANORM,
  final Box<double> RCOND,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
  final IPIV = IPIV_.having();
  final WORK = WORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  bool UPPER;
  int I, IP;
  final ISAVE = Array<int>(3);
  final AINVNM = Box(0.0);
  final KASE = Box(0);

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (ANORM < ZERO) {
    INFO.value = -5;
  }
  if (INFO.value != 0) {
    xerbla('ZHPCON', -INFO.value);
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

    IP = N * (N + 1) ~/ 2;
    for (I = N; I >= 1; I--) {
      if (IPIV[I] > 0 && AP[IP] == Complex.zero) return;
      IP -= I;
    }
  } else {
    // Lower triangular storage: examine D from top to bottom.

    IP = 1;
    for (I = 1; I <= N; I++) {
      if (IPIV[I] > 0 && AP[IP] == Complex.zero) return;
      IP += N - I + 1;
    }
  }

  // Estimate the 1-norm of the inverse.

  KASE.value = 0;
  while (true) {
    zlacn2(N, WORK(N + 1), WORK, AINVNM, KASE, ISAVE);
    if (KASE.value == 0) break;

    // Multiply by inv(L*D*L**H) or inv(U*D*U**H).

    zhptrs(UPLO, N, 1, AP, IPIV, WORK.asMatrix(), N, INFO);
  }

  // Compute the estimate of the reciprocal condition number.

  if (AINVNM.value != ZERO) RCOND.value = (ONE / AINVNM.value) / ANORM;
}
