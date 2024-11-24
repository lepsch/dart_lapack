// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zdrscl.dart';
import 'package:lapack/src/zlacn2.dart';
import 'package:lapack/src/zlantb.dart';
import 'package:lapack/src/zlatbs.dart';

void ztbcon(
  final String NORM,
  final String UPLO,
  final String DIAG,
  final int N,
  final int KD,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Box<double> RCOND,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
  final AB = AB_.having(ld: LDAB);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0, ZERO = 0.0;
  bool NOUNIT, ONENRM, UPPER;
  String NORMIN;
  int IX, KASE1;
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
  } else if (KD < 0) {
    INFO.value = -5;
  } else if (LDAB < KD + 1) {
    INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('ZTBCON', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) {
    RCOND.value = ONE;
    return;
  }

  RCOND.value = ZERO;
  SMLNUM = dlamch('Safe minimum') * max(N, 1);

  // Compute the 1-norm of the triangular matrix A or A**H.

  ANORM = zlantb(NORM, UPLO, DIAG, N, KD, AB, LDAB, RWORK);

  // Continue only if ANORM > 0.

  if (ANORM > ZERO) {
    // Estimate the 1-norm of the inverse of A.

    AINVNM.value = ZERO;
    NORMIN = 'N';
    if (ONENRM) {
      KASE1 = 1;
    } else {
      KASE1 = 2;
    }
    KASE.value = 0;
    while (true) {
      zlacn2(N, WORK(N + 1), WORK, AINVNM, KASE, ISAVE);
      if (KASE.value == 0) break;
      if (KASE.value == KASE1) {
        // Multiply by inv(A).

        zlatbs(UPLO, 'No transpose', DIAG, NORMIN, N, KD, AB, LDAB, WORK, SCALE,
            RWORK, INFO);
      } else {
        // Multiply by inv(A**H).

        zlatbs(UPLO, 'Conjugate transpose', DIAG, NORMIN, N, KD, AB, LDAB, WORK,
            SCALE, RWORK, INFO);
      }
      NORMIN = 'Y';

      // Multiply by 1/SCALE if doing so will not cause overflow.

      if (SCALE.value != ONE) {
        IX = izamax(N, WORK, 1);
        XNORM = WORK[IX].cabs1();
        if (SCALE.value < XNORM * SMLNUM || SCALE.value == ZERO) return;
        zdrscl(N, SCALE.value, WORK, 1);
      }
    }

    // Compute the estimate of the reciprocal condition number.

    if (AINVNM.value != ZERO) RCOND.value = (ONE / ANORM) / AINVNM.value;
  }
}
