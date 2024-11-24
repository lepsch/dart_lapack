// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/idamax.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacn2.dart';
import 'package:dart_lapack/src/dlantr.dart';
import 'package:dart_lapack/src/dlatrs.dart';
import 'package:dart_lapack/src/drscl.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dtrcon(
  final String NORM,
  final String UPLO,
  final String DIAG,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Box<double> RCOND,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ONE = 1.0, ZERO = 0.0;

  // Test the input parameters.

  INFO.value = 0;
  final UPPER = lsame(UPLO, 'U');
  final ONENRM = NORM == '1' || lsame(NORM, 'O');
  final NOUNIT = lsame(DIAG, 'N');

  if (!ONENRM && !lsame(NORM, 'I')) {
    INFO.value = -1;
  } else if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -2;
  } else if (!NOUNIT && !lsame(DIAG, 'U')) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (LDA < max(1, N)) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('DTRCON', -INFO.value);
    return;
  }

  // Quick return if possible
  if (N == 0) {
    RCOND.value = ONE;
    return;
  }

  RCOND.value = ZERO;
  final SMLNUM = dlamch('Safe minimum') * max(1, N);

  // Compute the norm of the triangular matrix A.
  final ANORM = dlantr(NORM, UPLO, DIAG, N, N, A, LDA, WORK);

  // Continue only if ANORM > 0.
  if (ANORM <= ZERO) return;

  // Estimate the norm of the inverse of A.
  final AINVNM = Box(ZERO);
  var NORMIN = 'N';
  final KASE1 = ONENRM ? 1 : 2;

  final ISAVE = Array<int>(3);
  final KASE = Box(0);
  while (true) {
    dlacn2(N, WORK(N + 1), WORK, IWORK, AINVNM, KASE, ISAVE);
    if (KASE.value == 0) break;

    final SCALE = Box(ZERO);
    if (KASE.value == KASE1) {
      // Multiply by inv(A).
      dlatrs(UPLO, 'No transpose', DIAG, NORMIN, N, A, LDA, WORK, SCALE,
          WORK(2 * N + 1), INFO);
    } else {
      // Multiply by inv(A**T).
      dlatrs(UPLO, 'Transpose', DIAG, NORMIN, N, A, LDA, WORK, SCALE,
          WORK(2 * N + 1), INFO);
    }
    NORMIN = 'Y';

    // Multiply by 1/SCALE if doing so will not cause overflow.
    if (SCALE.value != ONE) {
      final IX = idamax(N, WORK, 1);
      final XNORM = WORK[IX].abs();
      if (SCALE.value < XNORM * SMLNUM || SCALE.value == ZERO) return;
      drscl(N, SCALE.value, WORK, 1);
    }
  }

  // Compute the estimate of the reciprocal condition number.
  if (AINVNM.value != ZERO) RCOND.value = (ONE / ANORM) / AINVNM.value;
}
