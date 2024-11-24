// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlarzb.dart';
import 'package:lapack/src/dlarzt.dart';
import 'package:lapack/src/dlatrz.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dtzrzf(
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

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;
  bool LQUERY;
  int I,
      IB,
      IWS,
      KI,
      KK,
      LDWORK = 0,
      LWKMIN,
      LWKOPT = 0,
      M1,
      MU,
      NB = 0,
      NBMIN,
      NX;

  // Test the input arguments

  INFO.value = 0;
  LQUERY = (LWORK == -1);
  if (M < 0) {
    INFO.value = -1;
  } else if (N < M) {
    INFO.value = -2;
  } else if (LDA < max(1, M)) {
    INFO.value = -4;
  }

  if (INFO.value == 0) {
    if (M == 0 || M == N) {
      LWKOPT = 1;
      LWKMIN = 1;
    } else {
      // Determine the block size.

      NB = ilaenv(1, 'DGERQF', ' ', M, N, -1, -1);
      LWKOPT = M * NB;
      LWKMIN = max(1, M);
    }
    WORK[1] = LWKOPT.toDouble();

    if (LWORK < LWKMIN && !LQUERY) {
      INFO.value = -7;
    }
  }

  if (INFO.value != 0) {
    xerbla('DTZRZF', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (M == 0) {
    return;
  } else if (M == N) {
    for (I = 1; I <= N; I++) {
      TAU[I] = ZERO;
    }
    return;
  }

  NBMIN = 2;
  NX = 1;
  IWS = M;
  if (NB > 1 && NB < M) {
    // Determine when to cross over from blocked to unblocked code.

    NX = max(0, ilaenv(3, 'DGERQF', ' ', M, N, -1, -1));
    if (NX < M) {
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

  if (NB >= NBMIN && NB < M && NX < M) {
    // Use blocked code initially.
    // The last kk rows are handled by the block method.

    M1 = min(M + 1, N);
    KI = ((M - NX - 1) ~/ NB) * NB;
    KK = min(M, KI + NB);

    for (I = M - KK + KI + 1; I >= M - KK + 1; I -= NB) {
      IB = min(M - I + 1, NB);

      // Compute the TZ factorization of the current block
      // A(i:i+ib-1,i:n)

      dlatrz(IB, N - I + 1, N - M, A(I, I), LDA, TAU(I), WORK);
      if (I > 1) {
        // Form the triangular factor of the block reflector
        // H = H(i+ib-1) . . . H(i+1) H(i)

        dlarzt('Backward', 'Rowwise', N - M, IB, A(I, M1), LDA, TAU(I),
            WORK.asMatrix(LDWORK), LDWORK);

        // Apply H to A(1:i-1,i:n) from the right

        dlarzb(
            'Right',
            'No transpose',
            'Backward',
            'Rowwise',
            I - 1,
            N - I + 1,
            IB,
            N - M,
            A(I, M1),
            LDA,
            WORK.asMatrix(LDWORK),
            LDWORK,
            A(1, I),
            LDA,
            WORK(IB + 1).asMatrix(LDWORK),
            LDWORK);
      }
    }
    MU = I + NB - 1;
  } else {
    MU = M;
  }

  // Use unblocked code to factor the last or only block

  if (MU > 0) dlatrz(MU, N, N - M, A, LDA, TAU, WORK);

  WORK[1] = LWKOPT.toDouble();
}
