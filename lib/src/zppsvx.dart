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
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlanhp.dart';
import 'package:dart_lapack/src/zlaqhp.dart';
import 'package:dart_lapack/src/zppcon.dart';
import 'package:dart_lapack/src/zppequ.dart';
import 'package:dart_lapack/src/zpprfs.dart';
import 'package:dart_lapack/src/zpptrf.dart';
import 'package:dart_lapack/src/zpptrs.dart';

void zppsvx(
  final String FACT,
  final String UPLO,
  final int N,
  final int NRHS,
  final Array<Complex> AP_,
  final Array<Complex> AFP_,
  final Box<String> EQUED,
  final Array<double> S_,
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
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final AP = AP_.having();
  final AFP = AFP_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final S = S_.having();
  final FERR = FERR_.having();
  final BERR = BERR_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool EQUIL, NOFACT, RCEQU;
  int I, J;
  double ANORM, BIGNUM = 0, SMAX, SMIN, SMLNUM = 0;
  final AMAX = Box(0.0), SCOND = Box(0.0);
  final INFEQU = Box(0);

  INFO.value = 0;
  NOFACT = lsame(FACT, 'N');
  EQUIL = lsame(FACT, 'E');
  if (NOFACT || EQUIL) {
    EQUED.value = 'N';
    RCEQU = false;
  } else {
    RCEQU = lsame(EQUED.value, 'Y');
    SMLNUM = dlamch('Safe minimum');
    BIGNUM = ONE / SMLNUM;
  }

  // Test the input parameters.

  if (!NOFACT && !EQUIL && !lsame(FACT, 'F')) {
    INFO.value = -1;
  } else if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (NRHS < 0) {
    INFO.value = -4;
  } else if (lsame(FACT, 'F') && !(RCEQU || lsame(EQUED.value, 'N'))) {
    INFO.value = -7;
  } else {
    if (RCEQU) {
      SMIN = BIGNUM;
      SMAX = ZERO;
      for (J = 1; J <= N; J++) {
        SMIN = min(SMIN, S[J]);
        SMAX = max(SMAX, S[J]);
      }
      if (SMIN <= ZERO) {
        INFO.value = -8;
      } else if (N > 0) {
        SCOND.value = max(SMIN, SMLNUM) / min(SMAX, BIGNUM);
      } else {
        SCOND.value = ONE;
      }
    }
    if (INFO.value == 0) {
      if (LDB < max(1, N)) {
        INFO.value = -10;
      } else if (LDX < max(1, N)) {
        INFO.value = -12;
      }
    }
  }

  if (INFO.value != 0) {
    xerbla('ZPPSVX', -INFO.value);
    return;
  }

  if (EQUIL) {
    // Compute row and column scalings to equilibrate the matrix A.

    zppequ(UPLO, N, AP, S, SCOND, AMAX, INFEQU);
    if (INFEQU.value == 0) {
      // Equilibrate the matrix.

      zlaqhp(UPLO, N, AP, S, SCOND.value, AMAX.value, EQUED);
      RCEQU = lsame(EQUED.value, 'Y');
    }
  }

  // Scale the right-hand side.

  if (RCEQU) {
    for (J = 1; J <= NRHS; J++) {
      for (I = 1; I <= N; I++) {
        B[I][J] = S[I].toComplex() * B[I][J];
      }
    }
  }

  if (NOFACT || EQUIL) {
    // Compute the Cholesky factorization A = U**H * U or A = L * L**H.

    zcopy(N * (N + 1) ~/ 2, AP, 1, AFP, 1);
    zpptrf(UPLO, N, AFP, INFO);

    // Return if INFO is non-zero.

    if (INFO.value > 0) {
      RCOND.value = ZERO;
      return;
    }
  }

  // Compute the norm of the matrix A.

  ANORM = zlanhp('I', UPLO, N, AP, RWORK);

  // Compute the reciprocal of the condition number of A.

  zppcon(UPLO, N, AFP, ANORM, RCOND, WORK, RWORK, INFO);

  // Compute the solution matrix X.

  zlacpy('Full', N, NRHS, B, LDB, X, LDX);
  zpptrs(UPLO, N, NRHS, AFP, X, LDX, INFO);

  // Use iterative refinement to improve the computed solution and
  // compute error bounds and backward error estimates for it.

  zpprfs(UPLO, N, NRHS, AP, AFP, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO);

  // Transform the solution matrix X to a solution of the original
  // system.

  if (RCEQU) {
    for (J = 1; J <= NRHS; J++) {
      for (I = 1; I <= N; I++) {
        X[I][J] = S[I].toComplex() * X[I][J];
      }
    }
    for (J = 1; J <= NRHS; J++) {
      FERR[J] /= SCOND.value;
    }
  }

  // Set INFO = N+1 if the matrix is singular to working precision.

  if (RCOND.value < dlamch('Epsilon')) INFO.value = N + 1;
}
