// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlarzb.dart';
import 'package:dart_lapack/src/dlarzt.dart';
import 'package:dart_lapack/src/dormr3.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dormrz(
  final String SIDE,
  final String TRANS,
  final int M,
  final int N,
  final int K,
  final int L,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> TAU_,
  final Matrix<double> C_,
  final int LDC,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final C = C_.having(ld: LDC);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  const NBMAX = 64, LDT = NBMAX + 1, TSIZE = LDT * NBMAX;
  bool LEFT, LQUERY, NOTRAN;
  String TRANST;
  int I,
      I1,
      I2,
      I3,
      IB,
      IC = 0,
      IWT = 0,
      JA = 0,
      JC = 0,
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
  } else if (!NOTRAN && !lsame(TRANS, 'T')) {
    INFO.value = -2;
  } else if (M < 0) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (K < 0 || K > NQ) {
    INFO.value = -5;
  } else if (L < 0 || (LEFT && (L > M)) || (!LEFT && (L > N))) {
    INFO.value = -6;
  } else if (LDA < max(1, K)) {
    INFO.value = -8;
  } else if (LDC < max(1, M)) {
    INFO.value = -11;
  } else if (LWORK < NW && !LQUERY) {
    INFO.value = -13;
  }

  if (INFO.value == 0) {
    // Compute the workspace requirements

    if (M == 0 || N == 0) {
      LWKOPT = 1;
    } else {
      NB = min(NBMAX, ilaenv(1, 'DORMRQ', SIDE + TRANS, M, N, K, -1));
      LWKOPT = NW * NB + TSIZE;
    }
    WORK[1] = LWKOPT.toDouble();
  }

  if (INFO.value != 0) {
    xerbla('DORMRZ', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) {
    WORK[1] = 1;
    return;
  }

  NBMIN = 2;
  LDWORK = NW;
  if (NB > 1 && NB < K) {
    if (LWORK < LWKOPT) {
      NB = (LWORK - TSIZE) ~/ LDWORK;
      NBMIN = max(2, ilaenv(2, 'DORMRQ', SIDE + TRANS, M, N, K, -1));
    }
  }

  if (NB < NBMIN || NB >= K) {
    // Use unblocked code

    dormr3(SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC, WORK, IINFO);
  } else {
    // Use blocked code

    IWT = 1 + NW * NB;
    if ((LEFT && !NOTRAN) || (!LEFT && NOTRAN)) {
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
      JC = 1;
      JA = M - L + 1;
    } else {
      MI = M;
      IC = 1;
      JA = N - L + 1;
    }

    if (NOTRAN) {
      TRANST = 'T';
    } else {
      TRANST = 'N';
    }

    for (I = I1; I3 < 0 ? I >= I2 : I <= I2; I += I3) {
      IB = min(NB, K - I + 1);

      // Form the triangular factor of the block reflector
      // H = H(i+ib-1) . . . H(i+1) H(i)

      dlarzt('Backward', 'Rowwise', L, IB, A(I, JA), LDA, TAU(I),
          WORK(IWT).asMatrix(LDT), LDT);

      if (LEFT) {
        // H or H**T is applied to C(i:m,1:n)

        MI = M - I + 1;
        IC = I;
      } else {
        // H or H**T is applied to C(1:m,i:n)

        NI = N - I + 1;
        JC = I;
      }

      // Apply H or H**T

      dlarzb(
          SIDE,
          TRANST,
          'Backward',
          'Rowwise',
          MI,
          NI,
          IB,
          L,
          A(I, JA),
          LDA,
          WORK(IWT).asMatrix(LDT),
          LDT,
          C(IC, JC),
          LDC,
          WORK.asMatrix(LDWORK),
          LDWORK);
    }
  }

  WORK[1] = LWKOPT.toDouble();
}
