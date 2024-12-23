// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlaed1.dart';
import 'package:dart_lapack/src/dlaed7.dart';
import 'package:dart_lapack/src/dsteqr.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dlaed0(
  final int ICOMPQ,
  final int QSIZ,
  final int N,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<double> Q_,
  final int LDQ,
  final Matrix<double> QSTORE_,
  final int LDQS,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
  final D = D_.having();
  final E = E_.having();
  final Q = Q_.having(ld: LDQ);
  final QSTORE = QSTORE_.having(ld: LDQS);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  int CURLVL,
      CURPRB = 0,
      CURR,
      I,
      IGIVCL = 0,
      IGIVNM = 0,
      IGIVPT = 0,
      INDXQ,
      IPERM = 0,
      IPRMPT = 0,
      IQ = 0,
      IQPTR = 0,
      IWREM = 0,
      J,
      K,
      LGN,
      MATSIZ = 0,
      MSD2,
      SMLSIZ,
      SMM1,
      SPM1,
      SPM2,
      SUBMAT = 0,
      SUBPBS,
      TLVLS;
  double TEMP;

  // Test the input parameters.

  INFO.value = 0;

  if (ICOMPQ < 0 || ICOMPQ > 2) {
    INFO.value = -1;
  } else if ((ICOMPQ == 1) && (QSIZ < max(0, N))) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (LDQ < max(1, N)) {
    INFO.value = -7;
  } else if (LDQS < max(1, N)) {
    INFO.value = -9;
  }
  if (INFO.value != 0) {
    xerbla('DLAED0', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  SMLSIZ = ilaenv(9, 'DLAED0', ' ', 0, 0, 0, 0);

  // Determine the size and placement of the submatrices, and save in
  // the leading elements of IWORK.

  IWORK[1] = N;
  SUBPBS = 1;
  TLVLS = 0;
  while (IWORK[SUBPBS] > SMLSIZ) {
    for (J = SUBPBS; J >= 1; J--) {
      IWORK[2 * J] = (IWORK[J] + 1) ~/ 2;
      IWORK[2 * J - 1] = IWORK[J] ~/ 2;
    }
    TLVLS++;
    SUBPBS = 2 * SUBPBS;
  }
  for (J = 2; J <= SUBPBS; J++) {
    IWORK[J] += IWORK[J - 1];
  }

  // Divide the matrix into SUBPBS submatrices of size at most SMLSIZ+1
  // using rank-1 modifications (cuts).

  SPM1 = SUBPBS - 1;
  for (I = 1; I <= SPM1; I++) {
    SUBMAT = IWORK[I] + 1;
    SMM1 = SUBMAT - 1;
    D[SMM1] -= E[SMM1].abs();
    D[SUBMAT] -= E[SMM1].abs();
  }

  INDXQ = 4 * N + 3;
  if (ICOMPQ != 2) {
    // Set up workspaces for eigenvalues only/accumulate new vectors
    // routine

    TEMP = log(N) / log(TWO);
    LGN = TEMP.toInt();
    if (pow(2, LGN) < N) LGN++;
    if (pow(2, LGN) < N) LGN++;
    IPRMPT = INDXQ + N + 1;
    IPERM = IPRMPT + N * LGN;
    IQPTR = IPERM + N * LGN;
    IGIVPT = IQPTR + N + 2;
    IGIVCL = IGIVPT + N * LGN;

    IGIVNM = 1;
    IQ = IGIVNM + 2 * N * LGN;
    IWREM = IQ + pow(N, 2).toInt() + 1;

    // Initialize pointers

    for (I = 0; I <= SUBPBS; I++) {
      IWORK[IPRMPT + I] = 1;
      IWORK[IGIVPT + I] = 1;
    }
    IWORK[IQPTR] = 1;
  }

  // Solve each submatrix eigenproblem at the bottom of the divide and
  // conquer tree.

  CURR = 0;
  for (I = 0; I <= SPM1; I++) {
    if (I == 0) {
      SUBMAT = 1;
      MATSIZ = IWORK[1];
    } else {
      SUBMAT = IWORK[I] + 1;
      MATSIZ = IWORK[I + 1] - IWORK[I];
    }
    if (ICOMPQ == 2) {
      dsteqr('I', MATSIZ, D(SUBMAT), E(SUBMAT), Q(SUBMAT, SUBMAT), LDQ, WORK,
          INFO);
      if (INFO.value != 0) {
        INFO.value = SUBMAT * (N + 1) + SUBMAT + MATSIZ - 1;
        return;
      }
    } else {
      dsteqr(
          'I',
          MATSIZ,
          D(SUBMAT),
          E(SUBMAT),
          WORK(IQ - 1 + IWORK[IQPTR + CURR]).asMatrix(MATSIZ),
          MATSIZ,
          WORK,
          INFO);
      if (INFO.value != 0) {
        INFO.value = SUBMAT * (N + 1) + SUBMAT + MATSIZ - 1;
        return;
      }
      if (ICOMPQ == 1) {
        dgemm(
            'N',
            'N',
            QSIZ,
            MATSIZ,
            MATSIZ,
            ONE,
            Q(1, SUBMAT),
            LDQ,
            WORK(IQ - 1 + IWORK[IQPTR + CURR]).asMatrix(MATSIZ),
            MATSIZ,
            ZERO,
            QSTORE(1, SUBMAT),
            LDQS);
      }
      IWORK[IQPTR + CURR + 1] = IWORK[IQPTR + CURR] + pow(MATSIZ, 2).toInt();
      CURR++;
    }
    K = 1;
    for (J = SUBMAT; J <= IWORK[I + 1]; J++) {
      IWORK[INDXQ + J] = K;
      K++;
    }
  }

  // Successively merge eigensystems of adjacent submatrices
  // into eigensystem for the corresponding larger matrix.

  // while ( SUBPBS > 1 )

  CURLVL = 1;
  while (SUBPBS > 1) {
    SPM2 = SUBPBS - 2;
    for (I = 0; I <= SPM2; I += 2) {
      if (I == 0) {
        SUBMAT = 1;
        MATSIZ = IWORK[2];
        MSD2 = IWORK[1];
        CURPRB = 0;
      } else {
        SUBMAT = IWORK[I] + 1;
        MATSIZ = IWORK[I + 2] - IWORK[I];
        MSD2 = MATSIZ ~/ 2;
        CURPRB++;
      }

      // Merge lower order eigensystems (of size MSD2 and MATSIZ - MSD2)
      // into an eigensystem of size MATSIZ.
      // DLAED1 is used only for the full eigensystem of a tridiagonal
      // matrix.
      // DLAED7 handles the cases in which eigenvalues only or eigenvalues
      // and eigenvectors of a full symmetric matrix (which was reduced to
      // tridiagonal form) are desired.

      if (ICOMPQ == 2) {
        dlaed1(MATSIZ, D(SUBMAT), Q(SUBMAT, SUBMAT), LDQ, IWORK(INDXQ + SUBMAT),
            E.box(SUBMAT + MSD2 - 1), MSD2, WORK, IWORK(SUBPBS + 1), INFO);
      } else {
        dlaed7(
            ICOMPQ,
            MATSIZ,
            QSIZ,
            TLVLS,
            CURLVL,
            CURPRB,
            D(SUBMAT),
            QSTORE(1, SUBMAT),
            LDQS,
            IWORK(INDXQ + SUBMAT),
            E.box(SUBMAT + MSD2 - 1),
            MSD2,
            WORK(IQ),
            IWORK(IQPTR),
            IWORK(IPRMPT),
            IWORK(IPERM),
            IWORK(IGIVPT),
            IWORK(IGIVCL).asMatrix(2),
            WORK(IGIVNM).asMatrix(2),
            WORK(IWREM),
            IWORK(SUBPBS + 1),
            INFO);
      }
      if (INFO.value != 0) {
        INFO.value = SUBMAT * (N + 1) + SUBMAT + MATSIZ - 1;
        return;
      }
      IWORK[I ~/ 2 + 1] = IWORK[I + 2];
    }
    SUBPBS = SUBPBS ~/ 2;
    CURLVL++;
  }

  // Re-merge the eigenvalues/vectors which were deflated at the final
  // merge step.

  if (ICOMPQ == 1) {
    for (I = 1; I <= N; I++) {
      J = IWORK[INDXQ + I];
      WORK[I] = D[J];
      dcopy(QSIZ, QSTORE(1, J).asArray(), 1, Q(1, I).asArray(), 1);
    }
    dcopy(N, WORK, 1, D, 1);
  } else if (ICOMPQ == 2) {
    for (I = 1; I <= N; I++) {
      J = IWORK[INDXQ + I];
      WORK[I] = D[J];
      dcopy(N, Q(1, J).asArray(), 1, WORK(N * I + 1), 1);
    }
    dcopy(N, WORK, 1, D, 1);
    dlacpy('A', N, N, WORK(N + 1).asMatrix(N), N, Q, LDQ);
  } else {
    for (I = 1; I <= N; I++) {
      J = IWORK[INDXQ + I];
      WORK[I] = D[J];
    }
    dcopy(N, WORK, 1, D, 1);
  }
}
