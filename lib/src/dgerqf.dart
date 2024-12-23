// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgerq2.dart';
import 'package:dart_lapack/src/dlarfb.dart';
import 'package:dart_lapack/src/dlarft.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dgerqf(
  final int M,
  final int N,
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
  bool LQUERY;
  int I, IB, IWS, K = 0, KI, KK, LDWORK = 0, LWKOPT, MU, NB = 0, NBMIN, NU, NX;
  final IINFO = Box(0);

  // Test the input arguments

  INFO.value = 0;
  LQUERY = (LWORK == -1);
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, M)) {
    INFO.value = -4;
  }

  if (INFO.value == 0) {
    K = min(M, N);
    if (K == 0) {
      LWKOPT = 1;
    } else {
      NB = ilaenv(1, 'DGERQF', ' ', M, N, -1, -1);
      LWKOPT = M * NB;
    }
    WORK[1] = LWKOPT.toDouble();

    if (!LQUERY) {
      if (LWORK <= 0 || (N > 0 && LWORK < max(1, M))) INFO.value = -7;
    }
  }

  if (INFO.value != 0) {
    xerbla('DGERQF', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (K == 0) {
    return;
  }

  NBMIN = 2;
  NX = 1;
  IWS = M;
  if (NB > 1 && NB < K) {
    // Determine when to cross over from blocked to unblocked code.

    NX = max(0, ilaenv(3, 'DGERQF', ' ', M, N, -1, -1));
    if (NX < K) {
      // Determine if workspace is large enough for blocked code.

      LDWORK = M;
      IWS = LDWORK * NB;
      if (LWORK < IWS) {
        // Not enough workspace to use optimal NB:  reduce NB and
        // determine the minimum value of NB.

        NB = LWORK ~/ LDWORK;
        NBMIN = max(2, ilaenv(2, 'DGERQF', ' ', M, N, -1, -1));
      }
    }
  }

  if (NB >= NBMIN && NB < K && NX < K) {
    // Use blocked code initially.
    // The last kk rows are handled by the block method.

    KI = ((K - NX - 1) ~/ NB) * NB;
    KK = min(K, KI + NB);

    for (I = K - KK + KI + 1; I >= K - KK + 1; I -= NB) {
      IB = min(K - I + 1, NB);

      // Compute the RQ factorization of the current block
      // A(m-k+i:m-k+i+ib-1,1:n-k+i+ib-1)

      dgerq2(IB, N - K + I + IB - 1, A(M - K + I, 1), LDA, TAU(I), WORK, IINFO);
      if (M - K + I > 1) {
        // Form the triangular factor of the block reflector
        // H = H(i+ib-1) . . . H(i+1) H(i)

        dlarft('Backward', 'Rowwise', N - K + I + IB - 1, IB, A(M - K + I, 1),
            LDA, TAU(I), WORK.asMatrix(LDWORK), LDWORK);

        // Apply H to A(1:m-k+i-1,1:n-k+i+ib-1) from the right

        dlarfb(
            'Right',
            'No transpose',
            'Backward',
            'Rowwise',
            M - K + I - 1,
            N - K + I + IB - 1,
            IB,
            A(M - K + I, 1),
            LDA,
            WORK.asMatrix(LDWORK),
            LDWORK,
            A,
            LDA,
            WORK(IB + 1).asMatrix(LDWORK),
            LDWORK);
      }
    }
    MU = M - K + I + NB - 1;
    NU = N - K + I + NB - 1;
  } else {
    MU = M;
    NU = N;
  }

  // Use unblocked code to factor the last or only block

  if (MU > 0 && NU > 0) dgerq2(MU, NU, A, LDA, TAU, WORK, IINFO);

  WORK[1] = IWS.toDouble();
}
