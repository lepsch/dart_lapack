// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgttrs.dart';
import 'package:dart_lapack/src/dlacn2.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dgtcon(
  final String NORM,
  final int N,
  final Array<double> DL_,
  final Array<double> D_,
  final Array<double> DU_,
  final Array<double> DU2_,
  final Array<int> IPIV_,
  final double ANORM,
  final Box<double> RCOND,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
  final DL = DL_.having();
  final D = D_.having();
  final DU = DU_.having();
  final DU2 = DU2_.having();
  final IPIV = IPIV_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  bool ONENRM;
  int I, KASE1;
  final KASE = Box(0);
  final AINVNM = Box(0.0);
  final ISAVE = Array<int>(3);

  // Test the input arguments.

  INFO.value = 0;
  ONENRM = NORM == '1' || lsame(NORM, 'O');
  if (!ONENRM && !lsame(NORM, 'I')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (ANORM < ZERO) {
    INFO.value = -8;
  }
  if (INFO.value != 0) {
    xerbla('DGTCON', -INFO.value);
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

  // Check that D(1:N) is non-zero.

  for (I = 1; I <= N; I++) {
    if (D[I] == ZERO) return;
  }

  AINVNM.value = ZERO;
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
      // Multiply by inv(U)*inv(L).

      dgttrs('No transpose', N, 1, DL, D, DU, DU2, IPIV, WORK.asMatrix(N), N,
          INFO);
    } else {
      // Multiply by inv(L**T)*inv(U**T).

      dgttrs(
          'Transpose', N, 1, DL, D, DU, DU2, IPIV, WORK.asMatrix(N), N, INFO);
    }
  }
  // Compute the estimate of the reciprocal condition number.

  if (AINVNM.value != ZERO) RCOND.value = (ONE / AINVNM.value) / ANORM;
}
