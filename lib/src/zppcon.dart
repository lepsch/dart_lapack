// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/izamax.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zdrscl.dart';
import 'package:dart_lapack/src/zlacn2.dart';
import 'package:dart_lapack/src/zlatps.dart';

void zppcon(
  final String UPLO,
  final int N,
  final Array<Complex> AP_,
  final double ANORM,
  final Box<double> RCOND,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
  final AP = AP_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  bool UPPER;
  String NORMIN;
  int IX;
  double SCALE, SMLNUM;
  final ISAVE = Array<int>(3);
  final KASE = Box(0);
  final AINVNM = Box(0.0), SCALEL = Box(0.0), SCALEU = Box(0.0);

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
    xerbla('ZPPCON', -INFO.value);
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
    zlacn2(N, WORK(N + 1), WORK, AINVNM, KASE, ISAVE);
    if (KASE.value == 0) break;
    if (UPPER) {
      // Multiply by inv(U**H).

      zlatps('Upper', 'Conjugate transpose', 'Non-unit', NORMIN, N, AP, WORK,
          SCALEL, RWORK, INFO);
      NORMIN = 'Y';

      // Multiply by inv(U).

      zlatps('Upper', 'No transpose', 'Non-unit', NORMIN, N, AP, WORK, SCALEU,
          RWORK, INFO);
    } else {
      // Multiply by inv(L).

      zlatps('Lower', 'No transpose', 'Non-unit', NORMIN, N, AP, WORK, SCALEL,
          RWORK, INFO);
      NORMIN = 'Y';

      // Multiply by inv(L**H).

      zlatps('Lower', 'Conjugate transpose', 'Non-unit', NORMIN, N, AP, WORK,
          SCALEU, RWORK, INFO);
    }

    // Multiply by 1/SCALE if doing so will not cause overflow.

    SCALE = SCALEL.value * SCALEU.value;
    if (SCALE != ONE) {
      IX = izamax(N, WORK, 1);
      if (SCALE < WORK[IX].cabs1() * SMLNUM || SCALE == ZERO) return;
      zdrscl(N, SCALE, WORK, 1);
    }
  }

  // Compute the estimate of the reciprocal condition number.

  if (AINVNM.value != ZERO) RCOND.value = (ONE / AINVNM.value) / ANORM;
}
