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
import 'package:dart_lapack/src/zlarfb.dart';
import 'package:dart_lapack/src/zlarft.dart';
import 'package:dart_lapack/src/zunm2l.dart';

void zunmql(
  final String SIDE,
  final String TRANS,
  final int M,
  final int N,
  final int K,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> TAU_,
  final Matrix<Complex> C_,
  final int LDC,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final C = C_.having(ld: LDC);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  const NBMAX = 64, LDT = NBMAX + 1, TSIZE = LDT * NBMAX;
  bool LEFT, LQUERY, NOTRAN;
  int I,
      I1,
      I2,
      I3,
      IB,
      IWT,
      LDWORK,
      LWKOPT = 0,
      MI = 0,
      NB = 0,
      NBMIN,
      NI = 0,
      NQ,
      NW;
  final IINFO = Box(0);

  // Test the input arguments

  INFO.value = 0;
  LEFT = lsame(SIDE, 'L');
  NOTRAN = lsame(TRANS, 'N');
  LQUERY = (LWORK == -1);

  // NQ is the order of Q and NW is the minimum dimension of WORK

  if (LEFT) {
    NQ = M;
    NW = max(1, N);
  } else {
    NQ = N;
    NW = max(1, M);
  }
  if (!LEFT && !lsame(SIDE, 'R')) {
    INFO.value = -1;
  } else if (!NOTRAN && !lsame(TRANS, 'C')) {
    INFO.value = -2;
  } else if (M < 0) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (K < 0 || K > NQ) {
    INFO.value = -5;
  } else if (LDA < max(1, NQ)) {
    INFO.value = -7;
  } else if (LDC < max(1, M)) {
    INFO.value = -10;
  } else if (LWORK < NW && !LQUERY) {
    INFO.value = -12;
  }

  if (INFO.value == 0) {
    // Compute the workspace requirements

    if (M == 0 || N == 0) {
      LWKOPT = 1;
    } else {
      NB = min(NBMAX, ilaenv(1, 'ZUNMQL', SIDE + TRANS, M, N, K, -1));
      LWKOPT = NW * NB + TSIZE;
    }
    WORK[1] = LWKOPT.toComplex();
  }

  if (INFO.value != 0) {
    xerbla('ZUNMQL', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) {
    return;
  }

  NBMIN = 2;
  LDWORK = NW;
  if (NB > 1 && NB < K) {
    if (LWORK < LWKOPT) {
      NB = (LWORK - TSIZE) ~/ LDWORK;
      NBMIN = max(2, ilaenv(2, 'ZUNMQL', SIDE + TRANS, M, N, K, -1));
    }
  }

  if (NB < NBMIN || NB >= K) {
    // Use unblocked code

    zunm2l(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, IINFO);
  } else {
    // Use blocked code

    IWT = 1 + NW * NB;
    if ((LEFT && NOTRAN) || (!LEFT && !NOTRAN)) {
      I1 = 1;
      I2 = K;
      I3 = NB;
    } else {
      I1 = ((K - 1) ~/ NB) * NB + 1;
      I2 = 1;
      I3 = -NB;
    }

    if (LEFT) {
      NI = N;
    } else {
      MI = M;
    }

    for (I = I1; I3 < 0 ? I >= I2 : I <= I2; I += I3) {
      IB = min(NB, K - I + 1);

      // Form the triangular factor of the block reflector
      // H = H(i+ib-1) . . . H(i+1) H(i)

      zlarft('Backward', 'Columnwise', NQ - K + I + IB - 1, IB, A(1, I), LDA,
          TAU(I), WORK(IWT).asMatrix(), LDT);
      if (LEFT) {
        // H or H**H is applied to C(1:m-k+i+ib-1,1:n)

        MI = M - K + I + IB - 1;
      } else {
        // H or H**H is applied to C(1:m,1:n-k+i+ib-1)

        NI = N - K + I + IB - 1;
      }

      // Apply H or H**H

      zlarfb(SIDE, TRANS, 'Backward', 'Columnwise', MI, NI, IB, A(1, I), LDA,
          WORK(IWT).asMatrix(), LDT, C, LDC, WORK.asMatrix(), LDWORK);
    }
  }
  WORK[1] = LWKOPT.toComplex();
}
