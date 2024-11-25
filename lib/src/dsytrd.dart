// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dsyr2k.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlatrd.dart';
import 'package:dart_lapack/src/dsytd2.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dsytrd(
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> D_,
  final Array<double> E_,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final D = D_.having();
  final E = E_.having();
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  const ONE = 1.0;
  bool LQUERY, UPPER;
  int I, IWS, J, KK, LDWORK = 0, LWKOPT = 0, NB = 0, NBMIN, NX;
  final IINFO = Box(0);

  // Test the input parameters

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  LQUERY = LWORK == -1;
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  } else if (LWORK < 1 && !LQUERY) {
    INFO.value = -9;
  }

  if (INFO.value == 0) {
    // Determine the block size.

    NB = ilaenv(1, 'DSYTRD', UPLO, N, -1, -1, -1);
    LWKOPT = max(1, N * NB);
    WORK[1] = LWKOPT.toDouble();
  }

  if (INFO.value != 0) {
    xerbla('DSYTRD', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) {
    WORK[1] = 1;
    return;
  }

  NX = N;
  IWS = 1;
  if (NB > 1 && NB < N) {
    // Determine when to cross over from blocked to unblocked code
    // (last block is always handled by unblocked code).

    NX = max(NB, ilaenv(3, 'DSYTRD', UPLO, N, -1, -1, -1));
    if (NX < N) {
      // Determine if workspace is large enough for blocked code.

      LDWORK = N;
      IWS = LDWORK * NB;
      if (LWORK < IWS) {
        // Not enough workspace to use optimal NB:  determine the
        // minimum value of NB, and reduce NB or force use of
        // unblocked code by setting NX = N.

        NB = max(LWORK ~/ LDWORK, 1);
        NBMIN = ilaenv(2, 'DSYTRD', UPLO, N, -1, -1, -1);
        if (NB < NBMIN) NX = N;
      }
    } else {
      NX = N;
    }
  } else {
    NB = 1;
  }

  if (UPPER) {
    // Reduce the upper triangle of A.
    // Columns 1:kk are handled by the unblocked method.

    KK = N - ((N - NX + NB - 1) ~/ NB) * NB;
    for (I = N - NB + 1; I >= KK + 1; I -= NB) {
      // Reduce columns i:i+nb-1 to tridiagonal form and form the
      // matrix W which is needed to update the unreduced part of
      // the matrix

      dlatrd(
          UPLO, I + NB - 1, NB, A, LDA, E, TAU, WORK.asMatrix(LDWORK), LDWORK);

      // Update the unreduced submatrix A[1:i-1][1:i-1], using an
      // update of the form:  A := A - V*W**T - W*V**T

      dsyr2k(UPLO, 'No transpose', I - 1, NB, -ONE, A(1, I), LDA,
          WORK.asMatrix(LDWORK), LDWORK, ONE, A, LDA);

      // Copy superdiagonal elements back into A, and diagonal
      // elements into D

      for (J = I; J <= I + NB - 1; J++) {
        A[J - 1][J] = E[J - 1];
        D[J] = A[J][J];
      }
    }

    // Use unblocked code to reduce the last or only block

    dsytd2(UPLO, KK, A, LDA, D, E, TAU, IINFO);
  } else {
    // Reduce the lower triangle of A

    for (I = 1; I <= N - NX; I += NB) {
      // Reduce columns i:i+nb-1 to tridiagonal form and form the
      // matrix W which is needed to update the unreduced part of
      // the matrix

      dlatrd(UPLO, N - I + 1, NB, A(I, I), LDA, E(I), TAU(I),
          WORK.asMatrix(LDWORK), LDWORK);

      // Update the unreduced submatrix A[i+ib:n][i+ib:n], using
      // an update of the form:  A := A - V*W**T - W*V**T

      dsyr2k(UPLO, 'No transpose', N - I - NB + 1, NB, -ONE, A(I + NB, I), LDA,
          WORK(NB + 1).asMatrix(LDWORK), LDWORK, ONE, A(I + NB, I + NB), LDA);

      // Copy subdiagonal elements back into A, and diagonal
      // elements into D

      for (J = I; J <= I + NB - 1; J++) {
        A[J + 1][J] = E[J];
        D[J] = A[J][J];
      }
    }

    // Use unblocked code to reduce the last or only block

    dsytd2(UPLO, N - I + 1, A(I, I), LDA, D(I), E(I), TAU(I), IINFO);
  }

  WORK[1] = LWKOPT.toDouble();
}
