// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlamtsqr.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dorgtsqr(
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
  int LDC = 0, LWORKOPT = 0, LC = 0, LW = 0, NBLOCAL = 0, J;
  final IINFO = Box(0);

  // Test the input parameters

  LQUERY = LWORK == -1;
  INFO.value = 0;
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
  } else {
    // Test the input LWORK for the dimension of the array WORK.
    // This workspace is used to store array C(LDC, N) and WORK(LWORK)
    // in the call to DLAMTSQR. See the documentation for DLAMTSQR.

    if (LWORK < 2 && !LQUERY) {
      INFO.value = -10;
    } else {
      // Set block size for column blocks

      NBLOCAL = min(NB, N);

      // LWORK = -1, then set the size for the array C(LDC,N)
      // in DLAMTSQR call and set the optimal size of the work array
      // WORK(LWORK) in DLAMTSQR call.

      LDC = M;
      LC = LDC * N;
      LW = N * NBLOCAL;

      LWORKOPT = LC + LW;

      if ((LWORK < max(1, LWORKOPT)) && !LQUERY) {
        INFO.value = -10;
      }
    }
  }

  // Handle error in the input parameters and return workspace query.

  if (INFO.value != 0) {
    xerbla('DORGTSQR', -INFO.value);
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

  // (1) Form explicitly the tall-skinny M-by-N left submatrix Q1_in
  // of M-by-M orthogonal matrix Q_in, which is implicitly stored in
  // the subdiagonal part of input array A and in the input array T.
  // Perform by the following operation using the routine DLAMTSQR.

  // Q1_in = Q_in * ( I ), where I is a N-by-N identity matrix,
  //                ( 0 )        0 is a (M-N)-by-N zero matrix.

  // (1a) Form M-by-N matrix in the array WORK(1:LDC*N) with ones
  // on the diagonal and zeros elsewhere.

  dlaset('F', M, N, ZERO, ONE, WORK.asMatrix(LDC), LDC);

  // (1b)  On input, WORK(1:LDC*N) stores ( I );
  //                                      ( 0 )

  // On output, WORK(1:LDC*N) stores Q1_in.

  dlamtsqr('L', 'N', M, N, N, MB, NBLOCAL, A, LDA, T, LDT, WORK.asMatrix(LDC),
      LDC, WORK(LC + 1), LW, IINFO);

  // (2) Copy the result from the part of the work array (1:M,1:N)
  // with the leading dimension LDC that starts at WORK(1) into
  // the output array A(1:M,1:N) column-by-column.

  for (J = 1; J <= N; J++) {
    dcopy(M, WORK((J - 1) * LDC + 1), 1, A(1, J).asArray(), 1);
  }

  WORK[1] = LWORKOPT.toDouble();
}
