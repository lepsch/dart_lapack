// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlaed8.dart';
import 'package:dart_lapack/src/dlaed9.dart';
import 'package:dart_lapack/src/dlaeda.dart';
import 'package:dart_lapack/src/dlamrg.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dlaed7(
  final int ICOMPQ,
  final int N,
  final int QSIZ,
  final int TLVLS,
  final int CURLVL,
  final int CURPBM,
  final Array<double> D_,
  final Matrix<double> Q_,
  final int LDQ,
  final Array<int> INDXQ_,
  final Box<double> RHO,
  final int CUTPNT,
  final Array<double> QSTORE_,
  final Array<int> QPTR_,
  final Array<int> PRMPTR_,
  final Array<int> PERM_,
  final Array<int> GIVPTR_,
  final Matrix<int> GIVCOL_,
  final Matrix<double> GIVNUM_,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
  final D = D_.having();
  final Q = Q_.having(ld: LDQ);
  final INDXQ = INDXQ_.having();
  final QSTORE = QSTORE_.having();
  final QPTR = QPTR_.having();
  final PRMPTR = PRMPTR_.having();
  final PERM = PERM_.having();
  final GIVPTR = GIVPTR_.having();
  final GIVCOL = GIVCOL_.having(ld: 2);
  final GIVNUM = GIVNUM_.having(ld: 2);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  int COLTYP,
      CURR,
      I,
      IDLMDA,
      INDX,
      INDXC,
      INDXP,
      IQ2,
      IS,
      IW,
      IZ,
      LDQ2,
      N1,
      N2,
      PTR;
  final K = Box(0);

  // Test the input parameters.

  INFO.value = 0;

  if (ICOMPQ < 0 || ICOMPQ > 1) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (ICOMPQ == 1 && QSIZ < N) {
    INFO.value = -3;
  } else if (LDQ < max(1, N)) {
    INFO.value = -9;
  } else if (min(1, N) > CUTPNT || N < CUTPNT) {
    INFO.value = -12;
  }
  if (INFO.value != 0) {
    xerbla('DLAED7', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // The following values are for bookkeeping purposes only.  They are
  // integer pointers which indicate the portion of the workspace
  // used by a particular array in DLAED8 and DLAED9.

  if (ICOMPQ == 1) {
    LDQ2 = QSIZ;
  } else {
    LDQ2 = N;
  }

  IZ = 1;
  IDLMDA = IZ + N;
  IW = IDLMDA + N;
  IQ2 = IW + N;
  IS = IQ2 + N * LDQ2;

  INDX = 1;
  INDXC = INDX + N;
  COLTYP = INDXC + N;
  INDXP = COLTYP + N;

  // Form the z-vector which consists of the last row of Q_1 and the
  // first row of Q_2.

  PTR = 1 + pow(2, TLVLS).toInt();
  for (I = 1; I <= CURLVL - 1; I++) {
    PTR += pow(2, TLVLS - I).toInt();
  }
  CURR = PTR + CURPBM;
  dlaeda(N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR, GIVCOL, GIVNUM, QSTORE,
      QPTR, WORK(IZ), WORK(IZ + N), INFO);

  // When solving the final problem, we no longer need the stored data,
  // so we will overwrite the data from this level onto the previously
  // used storage space.

  if (CURLVL == TLVLS) {
    QPTR[CURR] = 1;
    PRMPTR[CURR] = 1;
    GIVPTR[CURR] = 1;
  }

  // Sort and Deflate eigenvalues.

  dlaed8(
      ICOMPQ,
      K,
      N,
      QSIZ,
      D,
      Q,
      LDQ,
      INDXQ,
      RHO,
      CUTPNT,
      WORK(IZ),
      WORK(IDLMDA),
      WORK(IQ2).asMatrix(LDQ2),
      LDQ2,
      WORK(IW),
      PERM(PRMPTR[CURR]),
      GIVPTR.box(CURR + 1),
      GIVCOL(1, GIVPTR[CURR]),
      GIVNUM(1, GIVPTR[CURR]),
      IWORK(INDXP),
      IWORK(INDX),
      INFO);
  PRMPTR[CURR + 1] = PRMPTR[CURR] + N;
  GIVPTR[CURR + 1] += GIVPTR[CURR];

  // Solve Secular Equation.

  if (K.value != 0) {
    dlaed9(
        K.value,
        1,
        K.value,
        N,
        D,
        WORK(IS).asMatrix(K.value),
        K.value,
        RHO.value,
        WORK(IDLMDA),
        WORK(IW),
        QSTORE(QPTR[CURR]).asMatrix(K.value),
        K.value,
        INFO);
    if (INFO.value != 0) return;
    if (ICOMPQ == 1) {
      dgemm('N', 'N', QSIZ, K.value, K.value, ONE, WORK(IQ2).asMatrix(LDQ2),
          LDQ2, QSTORE(QPTR[CURR]).asMatrix(K.value), K.value, ZERO, Q, LDQ);
    }
    QPTR[CURR + 1] = QPTR[CURR] + pow(K.value, 2).toInt();

    // Prepare the INDXQ sorting permutation.

    N1 = K.value;
    N2 = N - K.value;
    dlamrg(N1, N2, D, 1, -1, INDXQ);
  } else {
    QPTR[CURR + 1] = QPTR[CURR];
    for (I = 1; I <= N; I++) {
      INDXQ[I] = I;
    }
  }
}
