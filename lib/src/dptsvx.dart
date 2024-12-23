// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlanst.dart';
import 'package:dart_lapack/src/dptcon.dart';
import 'package:dart_lapack/src/dptrfs.dart';
import 'package:dart_lapack/src/dpttrf.dart';
import 'package:dart_lapack/src/dpttrs.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dptsvx(
  final String FACT,
  final int N,
  final int NRHS,
  final Array<double> D_,
  final Array<double> E_,
  final Array<double> DF_,
  final Array<double> EF_,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> X_,
  final int LDX,
  final Box<double> RCOND,
  final Array<double> FERR_,
  final Array<double> BERR_,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final D = D_.having();
  final E = E_.having();
  final DF = DF_.having();
  final EF = EF_.having();
  final FERR = FERR_.having();
  final BERR = BERR_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0;
  bool NOFACT;
  double ANORM;

  // Test the input parameters.

  INFO.value = 0;
  NOFACT = lsame(FACT, 'N');
  if (!NOFACT && !lsame(FACT, 'F')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDB < max(1, N)) {
    INFO.value = -9;
  } else if (LDX < max(1, N)) {
    INFO.value = -11;
  }
  if (INFO.value != 0) {
    xerbla('DPTSVX', -INFO.value);
    return;
  }

  if (NOFACT) {
    // Compute the L*D*L**T (or U**T*D*U) factorization of A.

    dcopy(N, D, 1, DF, 1);
    if (N > 1) dcopy(N - 1, E, 1, EF, 1);
    dpttrf(N, DF, EF, INFO);

    // Return if INFO is non-zero.

    if (INFO.value > 0) {
      RCOND.value = ZERO;
      return;
    }
  }

  // Compute the norm of the matrix A.

  ANORM = dlanst('1', N, D, E);

  // Compute the reciprocal of the condition number of A.

  dptcon(N, DF, EF, ANORM, RCOND, WORK, INFO);

  // Compute the solution vectors X.

  dlacpy('Full', N, NRHS, B, LDB, X, LDX);
  dpttrs(N, NRHS, DF, EF, X, LDX, INFO);

  // Use iterative refinement to improve the computed solutions and
  // compute error bounds and backward error estimates for them.

  dptrfs(N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR, BERR, WORK, INFO);

  // Set INFO = N+1 if the matrix is singular to working precision.

  if (RCOND.value < dlamch('Epsilon')) INFO.value = N + 1;
}
