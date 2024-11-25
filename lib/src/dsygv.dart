// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dtrmm.dart';
import 'package:dart_lapack/src/blas/dtrsm.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dpotrf.dart';
import 'package:dart_lapack/src/dsyev.dart';
import 'package:dart_lapack/src/dsygst.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dsygv(
  final int ITYPE,
  final String JOBZ,
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Array<double> W_,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final W = W_.having();
  final WORK = WORK_.having();
  const ONE = 1.0;
  bool LQUERY, UPPER, WANTZ;
  String TRANS;
  int LWKMIN, LWKOPT = 0, NB, NEIG;

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  UPPER = lsame(UPLO, 'U');
  LQUERY = (LWORK == -1);

  INFO.value = 0;
  if (ITYPE < 1 || ITYPE > 3) {
    INFO.value = -1;
  } else if (!(WANTZ || lsame(JOBZ, 'N'))) {
    INFO.value = -2;
  } else if (!(UPPER || lsame(UPLO, 'L'))) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (LDA < max(1, N)) {
    INFO.value = -6;
  } else if (LDB < max(1, N)) {
    INFO.value = -8;
  }

  if (INFO.value == 0) {
    LWKMIN = max(1, 3 * N - 1);
    NB = ilaenv(1, 'DSYTRD', UPLO, N, -1, -1, -1);
    LWKOPT = max(LWKMIN, (NB + 2) * N);
    WORK[1] = LWKOPT.toDouble();

    if (LWORK < LWKMIN && !LQUERY) {
      INFO.value = -11;
    }
  }

  if (INFO.value != 0) {
    xerbla('DSYGV', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Form a Cholesky factorization of B.

  dpotrf(UPLO, N, B, LDB, INFO);
  if (INFO.value != 0) {
    INFO.value = N + INFO.value;
    return;
  }

  // Transform problem to standard eigenvalue problem and solve.

  dsygst(ITYPE, UPLO, N, A, LDA, B, LDB, INFO);
  dsyev(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO);

  if (WANTZ) {
    // Backtransform eigenvectors to the original problem.

    NEIG = N;
    if (INFO.value > 0) NEIG = INFO.value - 1;
    if (ITYPE == 1 || ITYPE == 2) {
      // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
      // backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y

      if (UPPER) {
        TRANS = 'N';
      } else {
        TRANS = 'T';
      }

      dtrsm('Left', UPLO, TRANS, 'Non-unit', N, NEIG, ONE, B, LDB, A, LDA);
    } else if (ITYPE == 3) {
      // For B*A*x=(lambda)*x;
      // backtransform eigenvectors: x = L*y or U**T*y

      if (UPPER) {
        TRANS = 'T';
      } else {
        TRANS = 'N';
      }

      dtrmm('Left', UPLO, TRANS, 'Non-unit', N, NEIG, ONE, B, LDB, A, LDA);
    }
  }

  WORK[1] = LWKOPT.toDouble();
}
