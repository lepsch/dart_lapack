// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/idamax.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacn2.dart';
import 'package:dart_lapack/src/dlatps.dart';
import 'package:dart_lapack/src/drscl.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dppcon(
  final String UPLO,
  final int N,
  final Array<double> AP_,
  final double ANORM,
  final Box<double> RCOND,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
  final AP = AP_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  bool UPPER;
  String NORMIN;
  int IX;
  double SCALE, SMLNUM;
  final ISAVE = Array<int>(3);
  final AINVNM = Box(0.0), SCALEL = Box(0.0), SCALEU = Box(0.0);
  final KASE = Box(0);

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (ANORM < ZERO) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('DPPCON', -INFO.value);
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

  SMLNUM = dlamch('Safe minimum');

  // Estimate the 1-norm of the inverse.

  KASE.value = 0;
  NORMIN = 'N';
  while (true) {
    dlacn2(N, WORK(N + 1), WORK, IWORK, AINVNM, KASE, ISAVE);
    if (KASE.value == 0) break;
    if (UPPER) {
      // Multiply by inv(U**T).

      dlatps('Upper', 'Transpose', 'Non-unit', NORMIN, N, AP, WORK, SCALEL,
          WORK(2 * N + 1), INFO);
      NORMIN = 'Y';

      // Multiply by inv(U).

      dlatps('Upper', 'No transpose', 'Non-unit', NORMIN, N, AP, WORK, SCALEU,
          WORK(2 * N + 1), INFO);
    } else {
      // Multiply by inv(L).

      dlatps('Lower', 'No transpose', 'Non-unit', NORMIN, N, AP, WORK, SCALEL,
          WORK(2 * N + 1), INFO);
      NORMIN = 'Y';

      // Multiply by inv(L**T).

      dlatps('Lower', 'Transpose', 'Non-unit', NORMIN, N, AP, WORK, SCALEU,
          WORK(2 * N + 1), INFO);
    }

    // Multiply by 1/SCALE if doing so will not cause overflow.

    SCALE = SCALEL.value * SCALEU.value;
    if (SCALE != ONE) {
      IX = idamax(N, WORK, 1);
      if (SCALE < WORK[IX].abs() * SMLNUM || SCALE == ZERO) return;
      drscl(N, SCALE, WORK, 1);
    }
  }

  // Compute the estimate of the reciprocal condition number.

  if (AINVNM.value != ZERO) RCOND.value = (ONE / AINVNM.value) / ANORM;
}
