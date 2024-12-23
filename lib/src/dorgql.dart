// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlarfb.dart';
import 'package:dart_lapack/src/dlarft.dart';
import 'package:dart_lapack/src/dorg2l.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dorgql(
  final int M,
  final int N,
  final int K,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0;
  bool LQUERY;
  int I, IB, IWS, J, KK, L, LDWORK = 0, LWKOPT = 0, NB = 0, NBMIN, NX;
  final IINFO = Box(0);

  // Test the input arguments

  INFO.value = 0;
  LQUERY = (LWORK == -1);
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0 || N > M) {
    INFO.value = -2;
  } else if (K < 0 || K > N) {
    INFO.value = -3;
  } else if (LDA < max(1, M)) {
    INFO.value = -5;
  }

  if (INFO.value == 0) {
    if (N == 0) {
      LWKOPT = 1;
    } else {
      NB = ilaenv(1, 'DORGQL', ' ', M, N, K, -1);
      LWKOPT = N * NB;
    }
    WORK[1] = LWKOPT.toDouble();

    if (LWORK < max(1, N) && !LQUERY) {
      INFO.value = -8;
    }
  }

  if (INFO.value != 0) {
    xerbla('DORGQL', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N <= 0) {
    return;
  }

  NBMIN = 2;
  NX = 0;
  IWS = N;
  if (NB > 1 && NB < K) {
    // Determine when to cross over from blocked to unblocked code.

    NX = max(0, ilaenv(3, 'DORGQL', ' ', M, N, K, -1));
    if (NX < K) {
      // Determine if workspace is large enough for blocked code.

      LDWORK = N;
      IWS = LDWORK * NB;
      if (LWORK < IWS) {
        // Not enough workspace to use optimal NB:  reduce NB and
        // determine the minimum value of NB.

        NB = LWORK ~/ LDWORK;
        NBMIN = max(2, ilaenv(2, 'DORGQL', ' ', M, N, K, -1));
      }
    }
  }

  if (NB >= NBMIN && NB < K && NX < K) {
    // Use blocked code after the first block.
    // The last kk columns are handled by the block method.

    KK = min(K, ((K - NX + NB - 1) ~/ NB) * NB);

    // Set A(m-kk+1:m,1:n-kk) to zero.

    for (J = 1; J <= N - KK; J++) {
      for (I = M - KK + 1; I <= M; I++) {
        A[I][J] = ZERO;
      }
    }
  } else {
    KK = 0;
  }

  // Use unblocked code for the first or only block.

  dorg2l(M - KK, N - KK, K - KK, A, LDA, TAU, WORK, IINFO);

  if (KK > 0) {
    // Use blocked code

    for (I = K - KK + 1; I <= K; I += NB) {
      IB = min(NB, K - I + 1);
      if (N - K + I > 1) {
        // Form the triangular factor of the block reflector
        // H = H(i+ib-1) . . . H(i+1) H(i)

        dlarft('Backward', 'Columnwise', M - K + I + IB - 1, IB,
            A(1, N - K + I), LDA, TAU(I), WORK.asMatrix(LDWORK), LDWORK);

        // Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left

        dlarfb(
            'Left',
            'No transpose',
            'Backward',
            'Columnwise',
            M - K + I + IB - 1,
            N - K + I - 1,
            IB,
            A(1, N - K + I),
            LDA,
            WORK.asMatrix(LDWORK),
            LDWORK,
            A,
            LDA,
            WORK(IB + 1).asMatrix(LDWORK),
            LDWORK);
      }

      // Apply H to rows 1:m-k+i+ib-1 of current block

      dorg2l(M - K + I + IB - 1, IB, IB, A(1, N - K + I), LDA, TAU(I), WORK,
          IINFO);

      // Set rows m-k+i+ib:m of current block to zero

      for (J = N - K + I; J <= N - K + I + IB - 1; J++) {
        for (L = M - K + I + IB; L <= M; L++) {
          A[L][J] = ZERO;
        }
      }
    }
  }

  WORK[1] = IWS.toDouble();
}
