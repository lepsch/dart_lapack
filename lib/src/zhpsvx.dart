// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zcopy.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zhpcon.dart';
import 'package:dart_lapack/src/zhprfs.dart';
import 'package:dart_lapack/src/zhptrf.dart';
import 'package:dart_lapack/src/zhptrs.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlanhp.dart';

void zhpsvx(
  final String FACT,
  final String UPLO,
  final int N,
  final int NRHS,
  final Array<Complex> AP_,
  final Array<Complex> AFP_,
  final Array<int> IPIV_,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> X_,
  final int LDX,
  final Box<double> RCOND,
  final Array<double> FERR_,
  final Array<double> BERR_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
  final IPIV = IPIV_.having();
  final AP = AP_.having();
  final AFP = AFP_.having();
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final FERR = FERR_.having();
  final BERR = BERR_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const ZERO = 0.0;
  bool NOFACT;
  double ANORM;

  // Test the input parameters.

  INFO.value = 0;
  NOFACT = lsame(FACT, 'N');
  if (!NOFACT && !lsame(FACT, 'F')) {
    INFO.value = -1;
  } else if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (NRHS < 0) {
    INFO.value = -4;
  } else if (LDB < max(1, N)) {
    INFO.value = -9;
  } else if (LDX < max(1, N)) {
    INFO.value = -11;
  }
  if (INFO.value != 0) {
    xerbla('ZHPSVX', -INFO.value);
    return;
  }

  if (NOFACT) {
    // Compute the factorization A = U*D*U**H or A = L*D*L**H.

    zcopy(N * (N + 1) ~/ 2, AP, 1, AFP, 1);
    zhptrf(UPLO, N, AFP, IPIV, INFO);

    // Return if INFO is non-zero.

    if (INFO.value > 0) {
      RCOND.value = ZERO;
      return;
    }
  }

  // Compute the norm of the matrix A.

  ANORM = zlanhp('I', UPLO, N, AP, RWORK);

  // Compute the reciprocal of the condition number of A.

  zhpcon(UPLO, N, AFP, IPIV, ANORM, RCOND, WORK, INFO);

  // Compute the solution vectors X.

  zlacpy('Full', N, NRHS, B, LDB, X, LDX);
  zhptrs(UPLO, N, NRHS, AFP, IPIV, X, LDX, INFO);

  // Use iterative refinement to improve the computed solutions and
  // compute error bounds and backward error estimates for them.

  zhprfs(UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK,
      INFO);

  // Set INFO = N+1 if the matrix is singular to working precision.

  if (RCOND.value < dlamch('Epsilon')) INFO.value = N + 1;
}
