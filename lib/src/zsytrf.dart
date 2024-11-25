// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlasyf.dart';
import 'package:dart_lapack/src/zsytf2.dart';

void zsytrf(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final WORK = WORK_.having();
  bool LQUERY, UPPER;
  int IWS, J, K = 0, LDWORK, LWKOPT = 0, NB = 0, NBMIN;
  final KB = Box(0), IINFO = Box(0);

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  LQUERY = (LWORK == -1);
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  } else if (LWORK < 1 && !LQUERY) {
    INFO.value = -7;
  }

  if (INFO.value == 0) {
    // Determine the block size

    NB = ilaenv(1, 'ZSYTRF', UPLO, N, -1, -1, -1);
    LWKOPT = max(1, N * NB);
    WORK[1] = LWKOPT.toComplex();
  }

  if (INFO.value != 0) {
    xerbla('ZSYTRF', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  NBMIN = 2;
  LDWORK = N;
  if (NB > 1 && NB < N) {
    IWS = LDWORK * NB;
    if (LWORK < IWS) {
      NB = max(LWORK ~/ LDWORK, 1);
      NBMIN = max(2, ilaenv(2, 'ZSYTRF', UPLO, N, -1, -1, -1));
    }
  } else {
    IWS = 1;
  }
  if (NB < NBMIN) NB = N;

  if (UPPER) {
    // Factorize A as U*D*U**T using the upper triangle of A

    // K is the main loop index, decreasing from N to 1 in steps of
    // KB, where KB is the number of columns factorized by ZLASYF;
    // KB is either NB or NB-1, or K for the last block

    K = N;
    while (K >= 1) {
      if (K > NB) {
        // Factorize columns k-kb+1:k of A and use blocked code to
        // update columns 1:k-kb

        zlasyf(UPLO, K, NB, KB, A, LDA, IPIV, WORK.asMatrix(N), N, IINFO);
      } else {
        // Use unblocked code to factorize columns 1:k of A

        zsytf2(UPLO, K, A, LDA, IPIV, IINFO);
        KB.value = K;
      }

      // Set INFO on the first occurrence of a zero pivot

      if (INFO.value == 0 && IINFO.value > 0) INFO.value = IINFO.value;

      // Decrease K and return to the start of the main loop

      K -= KB.value;
    }
  } else {
    // Factorize A as L*D*L**T using the lower triangle of A

    // K is the main loop index, increasing from 1 to N in steps of
    // KB, where KB is the number of columns factorized by ZLASYF;
    // KB is either NB or NB-1, or N-K+1 for the last block

    K = 1;
    while (K <= N) {
      if (K <= N - NB) {
        // Factorize columns k:k+kb-1 of A and use blocked code to
        // update columns k+kb:n

        zlasyf(UPLO, N - K + 1, NB, KB, A(K, K), LDA, IPIV(K), WORK.asMatrix(N),
            N, IINFO);
      } else {
        // Use unblocked code to factorize columns k:n of A

        zsytf2(UPLO, N - K + 1, A(K, K), LDA, IPIV(K), IINFO);
        KB.value = N - K + 1;
      }

      // Set INFO on the first occurrence of a zero pivot

      if (INFO.value == 0 && IINFO.value > 0) INFO.value = IINFO.value + K - 1;

      // Adjust IPIV

      for (J = K; J <= K + KB.value - 1; J++) {
        if (IPIV[J] > 0) {
          IPIV[J] += K - 1;
        } else {
          IPIV[J] -= K - 1;
        }
      }

      // Increase K and return to the start of the main loop

      K += KB.value;
    }
  }
  WORK[1] = LWKOPT.toComplex();
}
