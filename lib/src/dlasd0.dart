// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlasd1.dart';
import 'package:dart_lapack/src/dlasdq.dart';
import 'package:dart_lapack/src/dlasdt.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dlasd0(
  final int N,
  final int SQRE,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<double> U_,
  final int LDU,
  final Matrix<double> VT_,
  final int LDVT,
  final int SMLSIZ,
  final Array<int> IWORK_,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
  final D = D_.having();
  final E = E_.having();
  final U = U_.having(ld: LDU);
  final VT = VT_.having(ld: LDVT);
  final IWORK = IWORK_.having();
  final WORK = WORK_.having();
  int I,
      I1,
      IC,
      IDXQ,
      IDXQC,
      IM1,
      INODE,
      ITEMP,
      IWK,
      J,
      LF,
      LL,
      LVL,
      M,
      NCC,
      NDB1,
      NDIML,
      NDIMR,
      NL,
      NLF,
      NLP1,
      NR,
      NRF,
      NRP1,
      SQREI;
  final ALPHA = Box(0.0), BETA = Box(0.0);
  final ND = Box(0), NLVL = Box(0);

  // Test the input parameters.

  INFO.value = 0;

  if (N < 0) {
    INFO.value = -1;
  } else if ((SQRE < 0) || (SQRE > 1)) {
    INFO.value = -2;
  }

  M = N + SQRE;

  if (LDU < N) {
    INFO.value = -6;
  } else if (LDVT < M) {
    INFO.value = -8;
  } else if (SMLSIZ < 3) {
    INFO.value = -9;
  }
  if (INFO.value != 0) {
    xerbla('DLASD0', -INFO.value);
    return;
  }

  // If the input matrix is too small, call DLASDQ to find the SVD.

  if (N <= SMLSIZ) {
    dlasdq('U', SQRE, N, M, N, 0, D, E, VT, LDVT, U, LDU, U, LDU, WORK, INFO);
    return;
  }

  // Set up the computation tree.

  INODE = 1;
  NDIML = INODE + N;
  NDIMR = NDIML + N;
  IDXQ = NDIMR + N;
  IWK = IDXQ + N;
  dlasdt(N, NLVL, ND, IWORK(INODE), IWORK(NDIML), IWORK(NDIMR), SMLSIZ);

  // For the nodes on bottom level of the tree, solve
  // their subproblems by DLASDQ.

  NDB1 = (ND.value + 1) ~/ 2;
  NCC = 0;
  for (I = NDB1; I <= ND.value; I++) {
    // IC : center row of each node
    // NL : number of rows of left  subproblem
    // NR : number of rows of right subproblem
    // NLF: starting row of the left   subproblem
    // NRF: starting row of the right  subproblem

    I1 = I - 1;
    IC = IWORK[INODE + I1];
    NL = IWORK[NDIML + I1];
    NLP1 = NL + 1;
    NR = IWORK[NDIMR + I1];
    NRP1 = NR + 1;
    NLF = IC - NL;
    NRF = IC + 1;
    SQREI = 1;
    dlasdq('U', SQREI, NL, NLP1, NL, NCC, D(NLF), E(NLF), VT(NLF, NLF), LDVT,
        U(NLF, NLF), LDU, U(NLF, NLF), LDU, WORK, INFO);
    if (INFO.value != 0) {
      return;
    }
    ITEMP = IDXQ + NLF - 2;
    for (J = 1; J <= NL; J++) {
      IWORK[ITEMP + J] = J;
    }
    if (I == ND.value) {
      SQREI = SQRE;
    } else {
      SQREI = 1;
    }
    NRP1 = NR + SQREI;
    dlasdq('U', SQREI, NR, NRP1, NR, NCC, D(NRF), E(NRF), VT(NRF, NRF), LDVT,
        U(NRF, NRF), LDU, U(NRF, NRF), LDU, WORK, INFO);
    if (INFO.value != 0) {
      return;
    }
    ITEMP = IDXQ + IC;
    for (J = 1; J <= NR; J++) {
      IWORK[ITEMP + J - 1] = J;
    }
  }

  // Now conquer each subproblem bottom-up.

  for (LVL = NLVL.value; LVL >= 1; LVL--) {
    // Find the first node LF and last node LL on the
    // current level LVL.

    if (LVL == 1) {
      LF = 1;
      LL = 1;
    } else {
      LF = pow(2, LVL - 1).toInt();
      LL = 2 * LF - 1;
    }
    for (I = LF; I <= LL; I++) {
      IM1 = I - 1;
      IC = IWORK[INODE + IM1];
      NL = IWORK[NDIML + IM1];
      NR = IWORK[NDIMR + IM1];
      NLF = IC - NL;
      if ((SQRE == 0) && (I == LL)) {
        SQREI = SQRE;
      } else {
        SQREI = 1;
      }
      IDXQC = IDXQ + NLF - 1;
      ALPHA.value = D[IC];
      BETA.value = E[IC];
      dlasd1(NL, NR, SQREI, D(NLF), ALPHA, BETA, U(NLF, NLF), LDU, VT(NLF, NLF),
          LDVT, IWORK(IDXQC), IWORK(IWK), WORK, INFO);

      // Report the possible convergence failure.

      if (INFO.value != 0) {
        return;
      }
    }
  }
}
