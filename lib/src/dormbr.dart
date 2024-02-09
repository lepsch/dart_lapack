import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dormlq.dart';
import 'package:lapack/src/dormqr.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dormbr(
  final String VECT,
  final String SIDE,
  final String TRANS,
  final int M,
  final int N,
  final int K,
  final Matrix<double> A,
  final int LDA,
  final Array<double> TAU,
  final Matrix<double> C,
  final int LDC,
  final Array<double> WORK,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  bool APPLYQ, LEFT, LQUERY, NOTRAN;
  String TRANST;
  int I1, I2, LWKOPT = 0, MI, NB, NI, NQ, NW;
  final IINFO = Box(0);

  // Test the input arguments

  INFO.value = 0;
  APPLYQ = lsame(VECT, 'Q');
  LEFT = lsame(SIDE, 'L');
  NOTRAN = lsame(TRANS, 'N');
  LQUERY = (LWORK == -1);

  // NQ is the order of Q or P and NW is the minimum dimension of WORK

  if (LEFT) {
    NQ = M;
    NW = max(1, N);
  } else {
    NQ = N;
    NW = max(1, M);
  }
  if (!APPLYQ && !lsame(VECT, 'P')) {
    INFO.value = -1;
  } else if (!LEFT && !lsame(SIDE, 'R')) {
    INFO.value = -2;
  } else if (!NOTRAN && !lsame(TRANS, 'T')) {
    INFO.value = -3;
  } else if (M < 0) {
    INFO.value = -4;
  } else if (N < 0) {
    INFO.value = -5;
  } else if (K < 0) {
    INFO.value = -6;
  } else if ((APPLYQ && LDA < max(1, NQ)) ||
      (!APPLYQ && LDA < max(1, min(NQ, K)))) {
    INFO.value = -8;
  } else if (LDC < max(1, M)) {
    INFO.value = -11;
  } else if (LWORK < NW && !LQUERY) {
    INFO.value = -13;
  }

  if (INFO.value == 0) {
    if (APPLYQ) {
      if (LEFT) {
        NB = ilaenv(1, 'DORMQR', SIDE + TRANS, M - 1, N, M - 1, -1);
      } else {
        NB = ilaenv(1, 'DORMQR', SIDE + TRANS, M, N - 1, N - 1, -1);
      }
    } else {
      if (LEFT) {
        NB = ilaenv(1, 'DORMLQ', SIDE + TRANS, M - 1, N, M - 1, -1);
      } else {
        NB = ilaenv(1, 'DORMLQ', SIDE + TRANS, M, N - 1, N - 1, -1);
      }
    }
    LWKOPT = NW * NB;
    WORK[1] = LWKOPT.toDouble();
  }

  if (INFO.value != 0) {
    xerbla('DORMBR', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  WORK[1] = 1;
  if (M == 0 || N == 0) return;

  if (APPLYQ) {
    // Apply Q

    if (NQ >= K) {
      // Q was determined by a call to DGEBRD with nq >= k

      dormqr(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, IINFO);
    } else if (NQ > 1) {
      // Q was determined by a call to DGEBRD with nq < k

      if (LEFT) {
        MI = M - 1;
        NI = N;
        I1 = 2;
        I2 = 1;
      } else {
        MI = M;
        NI = N - 1;
        I1 = 1;
        I2 = 2;
      }
      dormqr(
        SIDE,
        TRANS,
        MI,
        NI,
        NQ - 1,
        A(2,1),
        LDA,
        TAU,
        C(I1,I2),
        LDC,
        WORK,
        LWORK,
        IINFO,
      );
    }
  } else {
    // Apply P

    if (NOTRAN) {
      TRANST = 'T';
    } else {
      TRANST = 'N';
    }
    if (NQ > K) {
      // P was determined by a call to DGEBRD with nq > k

      dormlq(SIDE, TRANST, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, IINFO);
    } else if (NQ > 1) {
      // P was determined by a call to DGEBRD with nq <= k

      if (LEFT) {
        MI = M - 1;
        NI = N;
        I1 = 2;
        I2 = 1;
      } else {
        MI = M;
        NI = N - 1;
        I1 = 1;
        I2 = 2;
      }
      dormlq(
        SIDE,
        TRANST,
        MI,
        NI,
        NQ - 1,
        A(1, 2),
        LDA,
        TAU,
        C(I1, I2),
        LDC,
        WORK,
        LWORK,
        IINFO,
      );
    }
  }
  WORK[1] = LWKOPT.toDouble();
  return;
}
