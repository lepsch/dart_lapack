// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlatsqr.dart';
import 'package:dart_lapack/src/dorgtsqr_row.dart';
import 'package:dart_lapack/src/dorhr_col.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dgetsqrhrt(
  final int M,
  final int N,
  final int MB1,
  final int NB1,
  final int NB2,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> T_,
  final int LDT,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final T = T_.having(ld: LDT);
  final WORK = WORK_.having();
  const ONE = 1.0;
  bool LQUERY;
  int I,
      J,
      LW1 = 0,
      LW2 = 0,
      LWT = 0,
      LDWT = 0,
      LWORKOPT = 0,
      NB1LOCAL = 0,
      NB2LOCAL,
      NUM_ALL_ROW_BLOCKS;
  final IINFO = Box(0);
  // Test the input arguments

  INFO.value = 0;
  LQUERY = (LWORK == -1);
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0 || M < N) {
    INFO.value = -2;
  } else if (MB1 <= N) {
    INFO.value = -3;
  } else if (NB1 < 1) {
    INFO.value = -4;
  } else if (NB2 < 1) {
    INFO.value = -5;
  } else if (LDA < max(1, M)) {
    INFO.value = -7;
  } else if (LDT < max(1, min(NB2, N))) {
    INFO.value = -9;
  } else {
    // Test the input LWORK for the dimension of the array WORK.
    // This workspace is used to store array:
    // a) Matrix T and WORK for DLATSQR;
    // b) N-by-N upper-triangular factor R_tsqr;
    // c) Matrix T and array WORK for DORGTSQR_ROW;
    // d) Diagonal D for DORHR_COL.

    if (LWORK < N * N + 1 && !LQUERY) {
      INFO.value = -11;
    } else {
      // Set block size for column blocks

      NB1LOCAL = min(NB1, N);

      NUM_ALL_ROW_BLOCKS = max(1, ((M - N) / (MB1 - N)).ceil());

      // Length and leading dimension of WORK array to place
      // T array in TSQR.

      LWT = NUM_ALL_ROW_BLOCKS * N * NB1LOCAL;

      LDWT = NB1LOCAL;

      // Length of TSQR work array

      LW1 = NB1LOCAL * N;

      // Length of DORGTSQR_ROW work array.

      LW2 = NB1LOCAL * max(NB1LOCAL, (N - NB1LOCAL));

      LWORKOPT = max(LWT + LW1, max(LWT + N * N + LW2, LWT + N * N + N));
      LWORKOPT = max(1, LWORKOPT);

      if (LWORK < LWORKOPT && !LQUERY) {
        INFO.value = -11;
      }
    }
  }

  // Handle error in the input parameters and return workspace query.

  if (INFO.value != 0) {
    xerbla('DGETSQRHRT', -INFO.value);
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

  NB2LOCAL = min(NB2, N);

  // (1) Perform TSQR-factorization of the M-by-N matrix A.

  dlatsqr(M, N, MB1, NB1LOCAL, A, LDA, WORK.asMatrix(LDWT), LDWT, WORK(LWT + 1),
      LW1, IINFO);

  // (2) Copy the factor R_tsqr stored in the upper-triangular part
  //     of A into the square matrix in the work array
  //     WORK(LWT+1:LWT+N*N) column-by-column.

  for (J = 1; J <= N; J++) {
    dcopy(J, A(1, J).asArray(), 1, WORK(LWT + N * (J - 1) + 1), 1);
  }

  // (3) Generate a M-by-N matrix Q with orthonormal columns from
  // the result stored below the diagonal in the array A in place.

  dorgtsqr_row(M, N, MB1, NB1LOCAL, A, LDA, WORK.asMatrix(LDWT), LDWT,
      WORK(LWT + N * N + 1), LW2, IINFO);

  // (4) Perform the reconstruction of Householder vectors from
  // the matrix Q (stored in A) in place.

  dorhr_col(M, N, NB2LOCAL, A, LDA, T, LDT, WORK(LWT + N * N + 1), IINFO);

  // (5) Copy the factor R_tsqr stored in the square matrix in the
  // work array WORK(LWT+1:LWT+N*N) into the upper-triangular
  // part of A.

  // (6) Compute from R_tsqr the factor R_hr corresponding to
  // the reconstructed Householder vectors, i.e. R_hr = S * R_tsqr.
  // This multiplication by the sign matrix S on the left means
  // changing the sign of I-th row of the matrix R_tsqr according
  // to sign of the I-th diagonal element DIAG(I) of the matrix S.
  // DIAG is stored in WORK( LWT+N*N+1 ) from the DORHR_COL output.

  // (5) and (6) can be combined in a single loop, so the rows in A
  // are accessed only once.

  for (I = 1; I <= N; I++) {
    if (WORK[LWT + N * N + I] == -ONE) {
      for (J = I; J <= N; J++) {
        A[I][J] = -ONE * WORK[LWT + N * (J - 1) + I];
      }
    } else {
      dcopy(N - I + 1, WORK(LWT + N * (I - 1) + I), N, A(I, I).asArray(), LDA);
    }
  }

  WORK[1] = LWORKOPT.toDouble();
}
