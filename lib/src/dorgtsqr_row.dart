// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlarfb_gett.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dorgtsqr_row(
  final int M,
  final int N,
  final int MB,
  final int NB,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> T_,
  final int LDT,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final T = T_.having(ld: LDT);
  final WORK = WORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  bool LQUERY;
  int NBLOCAL,
      MB2,
      M_PLUS_ONE,
      ITMP,
      IB_BOTTOM,
      LWORKOPT = 0,
      NUM_ALL_ROW_BLOCKS,
      JB_T,
      IB,
      IMB,
      KB,
      KB_LAST,
      KNB,
      MB1;
  final DUMMY = Matrix<double>(1, 1);

  // Test the input parameters

  INFO.value = 0;
  LQUERY = LWORK == -1;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0 || M < N) {
    INFO.value = -2;
  } else if (MB <= N) {
    INFO.value = -3;
  } else if (NB < 1) {
    INFO.value = -4;
  } else if (LDA < max(1, M)) {
    INFO.value = -6;
  } else if (LDT < max(1, min(NB, N))) {
    INFO.value = -8;
  } else if (LWORK < 1 && !LQUERY) {
    INFO.value = -10;
  }

  NBLOCAL = min(NB, N);

  // Determine the workspace size.

  if (INFO.value == 0) {
    LWORKOPT = NBLOCAL * max(NBLOCAL, (N - NBLOCAL));
  }

  // Handle error in the input parameters and handle the workspace query.

  if (INFO.value != 0) {
    xerbla('DORGTSQR_ROW', -INFO.value);
    return;
  } else if (LQUERY) {
    WORK[1] = LWORKOPT.toDouble();
    return;
  }

  // Quick return if possible

  if (min(M, N) == 0) {
    WORK[1] = LWORKOPT.toDouble();
    return;
  }

  // (0) Set the upper-triangular part of the matrix A to zero and
  // its diagonal elements to one.

  dlaset('U', M, N, ZERO, ONE, A, LDA);

  // KB_LAST is the column index of the last column block reflector
  // in the matrices T and V.

  KB_LAST = ((N - 1) ~/ NBLOCAL) * NBLOCAL + 1;

  // (1) Bottom-up loop over row blocks of A, except the top row block.
  // NOTE: If MB>=M, then the loop is never executed.

  if (MB < M) {
    // MB2 is the row blocking size for the row blocks before the
    // first top row block in the matrix A. IB is the row index for
    // the row blocks in the matrix A before the first top row block.
    // IB_BOTTOM is the row index for the last bottom row block
    // in the matrix A. JB_T is the column index of the corresponding
    // column block in the matrix T.

    // Initialize variables.

    // NUM_ALL_ROW_BLOCKS is the number of row blocks in the matrix A
    // including the first row block.

    MB2 = MB - N;
    M_PLUS_ONE = M + 1;
    ITMP = (M - MB - 1) ~/ MB2;
    IB_BOTTOM = ITMP * MB2 + MB + 1;
    NUM_ALL_ROW_BLOCKS = ITMP + 2;
    JB_T = NUM_ALL_ROW_BLOCKS * N + 1;

    for (IB = IB_BOTTOM; IB >= MB + 1; IB -= MB2) {
      // Determine the block size IMB for the current row block
      // in the matrix A.

      IMB = min(M_PLUS_ONE - IB, MB2);

      // Determine the column index JB_T for the current column block
      // in the matrix T.

      JB_T -= N;

      // Apply column blocks of H in the row block from right to left.

      // KB is the column index of the current column block reflector
      // in the matrices T and V.

      for (KB = KB_LAST; KB >= 1; KB -= NBLOCAL) {
        // Determine the size of the current column block KNB in
        // the matrices T and V.

        KNB = min(NBLOCAL, N - KB + 1);

        dlarfb_gett('I', IMB, N - KB + 1, KNB, T(1, JB_T + KB - 1), LDT,
            A(KB, KB), LDA, A(IB, KB), LDA, WORK.asMatrix(KNB), KNB);
      }
    }
  }

  // (2) Top row block of A.
  // NOTE: If MB>=M, then we have only one row block of A of size M
  // and we work on the entire matrix A.

  MB1 = min(MB, M);

  // Apply column blocks of H in the top row block from right to left.

  // KB is the column index of the current block reflector in
  // the matrices T and V.

  for (KB = KB_LAST; KB >= 1; KB -= NBLOCAL) {
    // Determine the size of the current column block KNB in
    // the matrices T and V.

    KNB = min(NBLOCAL, N - KB + 1);

    if (MB1 - KB - KNB + 1 == 0) {
      // In SLARFB_GETT parameters, when M=0, then the matrix B
      // does not exist, hence we need to pass a dummy array
      // reference DUMMY(1,1) to B with LDDUMMY=1.

      dlarfb_gett('N', 0, N - KB + 1, KNB, T(1, KB), LDT, A(KB, KB), LDA,
          DUMMY(1, 1), 1, WORK.asMatrix(KNB), KNB);
    } else {
      dlarfb_gett('N', MB1 - KB - KNB + 1, N - KB + 1, KNB, T(1, KB), LDT,
          A(KB, KB), LDA, A(KB + KNB, KB), LDA, WORK.asMatrix(KNB), KNB);
    }
  }

  WORK[1] = LWORKOPT.toDouble();
}
