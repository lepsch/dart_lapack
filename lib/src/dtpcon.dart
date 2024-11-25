// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/idamax.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacn2.dart';
import 'package:dart_lapack/src/dlantp.dart';
import 'package:dart_lapack/src/dlatps.dart';
import 'package:dart_lapack/src/drscl.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dtpcon(
  final String NORM,
  final String UPLO,
  final String DIAG,
  final int N,
  final Array<double> AP_,
  final Box<double> RCOND,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
  final AP = AP_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  bool NOUNIT, ONENRM, UPPER;
  String NORMIN;
  int IX, KASE1 = 0;
  double ANORM, SMLNUM, XNORM;
  final ISAVE = Array<int>(3);
  final AINVNM = Box(0.0), SCALE = Box(0.0);
  final KASE = Box(0);

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  ONENRM = NORM == '1' || lsame(NORM, 'O');
  NOUNIT = lsame(DIAG, 'N');

  if (!ONENRM && !lsame(NORM, 'I')) {
    INFO.value = -1;
  } else if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -2;
  } else if (!NOUNIT && !lsame(DIAG, 'U')) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('DTPCON', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) {
    RCOND.value = ONE;
    return;
  }

  RCOND.value = ZERO;
  SMLNUM = dlamch('Safe minimum') * max(1, N);

  // Compute the norm of the triangular matrix A.

  ANORM = dlantp(NORM, UPLO, DIAG, N, AP, WORK);

  // Continue only if ANORM > 0.

  if (ANORM > ZERO) {
    // Estimate the norm of the inverse of A.

    AINVNM.value = ZERO;
    NORMIN = 'N';
    if (ONENRM) {
      KASE1 = 1;
    } else {
      KASE1 = 2;
    }
    KASE.value = 0;
    while (true) {
      dlacn2(N, WORK(N + 1), WORK, IWORK, AINVNM, KASE, ISAVE);
      if (KASE.value == 0) break;
      if (KASE.value == KASE1) {
        // Multiply by inv(A).

        dlatps(UPLO, 'No transpose', DIAG, NORMIN, N, AP, WORK, SCALE,
            WORK(2 * N + 1), INFO);
      } else {
        // Multiply by inv(A**T).

        dlatps(UPLO, 'Transpose', DIAG, NORMIN, N, AP, WORK, SCALE,
            WORK(2 * N + 1), INFO);
      }
      NORMIN = 'Y';

      // Multiply by 1/SCALE if doing so will not cause overflow.

      if (SCALE.value != ONE) {
        IX = idamax(N, WORK, 1);
        XNORM = WORK[IX].abs();
        if (SCALE.value < XNORM * SMLNUM || SCALE.value == ZERO) return;
        drscl(N, SCALE.value, WORK, 1);
      }
    }

    // Compute the estimate of the reciprocal condition number.

    if (AINVNM.value != ZERO) RCOND.value = (ONE / ANORM) / AINVNM.value;
  }
}
