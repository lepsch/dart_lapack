// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlaed2.dart';
import 'package:dart_lapack/src/dlaed3.dart';
import 'package:dart_lapack/src/dlamrg.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dlaed1(
  final int N,
  final Array<double> D_,
  final Matrix<double> Q_,
  final int LDQ,
  final Array<int> INDXQ_,
  final Box<double> RHO,
  final int CUTPNT,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
  final D = D_.having();
  final Q = Q_.having(ld: LDQ);
  final INDXQ = INDXQ_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  int COLTYP, I, IDLMDA, INDX, INDXC, INDXP, IQ2, IS, IW, IZ, N1, N2, ZPP1;
  final K = Box(0);

  // Test the input parameters.

  INFO.value = 0;

  if (N < 0) {
    INFO.value = -1;
  } else if (LDQ < max(1, N)) {
    INFO.value = -4;
  } else if (min(1, N ~/ 2) > CUTPNT || (N ~/ 2) < CUTPNT) {
    INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('DLAED1', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // The following values are integer pointers which indicate
  // the portion of the workspace
  // used by a particular array in DLAED2 and DLAED3.

  IZ = 1;
  IDLMDA = IZ + N;
  IW = IDLMDA + N;
  IQ2 = IW + N;

  INDX = 1;
  INDXC = INDX + N;
  COLTYP = INDXC + N;
  INDXP = COLTYP + N;

  // Form the z-vector which consists of the last row of Q_1 and the
  // first row of Q_2.

  dcopy(CUTPNT, Q(CUTPNT, 1).asArray(), LDQ, WORK(IZ), 1);
  ZPP1 = CUTPNT + 1;
  dcopy(N - CUTPNT, Q(ZPP1, ZPP1).asArray(), LDQ, WORK(IZ + CUTPNT), 1);

  // Deflate eigenvalues.

  dlaed2(K, N, CUTPNT, D, Q, LDQ, INDXQ, RHO, WORK(IZ), WORK(IDLMDA), WORK(IW),
      WORK(IQ2), IWORK(INDX), IWORK(INDXC), IWORK(INDXP), IWORK(COLTYP), INFO);

  if (INFO.value != 0) return;

  // Solve Secular Equation.

  if (K.value != 0) {
    IS = (IWORK[COLTYP] + IWORK[COLTYP + 1]) * CUTPNT +
        (IWORK[COLTYP + 1] + IWORK[COLTYP + 2]) * (N - CUTPNT) +
        IQ2;
    dlaed3(K.value, N, CUTPNT, D, Q, LDQ, RHO.value, WORK(IDLMDA), WORK(IQ2),
        IWORK(INDXC), IWORK(COLTYP), WORK(IW), WORK(IS), INFO);
    if (INFO.value != 0) return;

    // Prepare the INDXQ sorting permutation.

    N1 = K.value;
    N2 = N - K.value;
    dlamrg(N1, N2, D, 1, -1, INDXQ);
  } else {
    for (I = 1; I <= N; I++) {
      INDXQ[I] = I;
    }
  }
}
