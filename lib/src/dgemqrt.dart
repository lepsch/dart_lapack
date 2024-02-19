import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlarfb.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgemqrt(
  final String SIDE,
  final String TRANS,
  final int M,
  final int N,
  final int K,
  final int NB,
  final Matrix<double> V_,
  final int LDV,
  final Matrix<double> T_,
  final int LDT,
  final Matrix<double> C_,
  final int LDC,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final V = V_.dim(LDV);
  final T = T_.dim(LDT);
  final C = C_.dim(LDC);
  final WORK = WORK_.dim();
  bool LEFT, RIGHT, TRAN, NOTRAN;
  int I, IB, LDWORK = 0, KF, Q = 0;

  // .. Test the input arguments ..

  INFO.value = 0;
  LEFT = lsame(SIDE, 'L');
  RIGHT = lsame(SIDE, 'R');
  TRAN = lsame(TRANS, 'T');
  NOTRAN = lsame(TRANS, 'N');

  if (LEFT) {
    LDWORK = max(1, N);
    Q = M;
  } else if (RIGHT) {
    LDWORK = max(1, M);
    Q = N;
  }
  if (!LEFT && !RIGHT) {
    INFO.value = -1;
  } else if (!TRAN && !NOTRAN) {
    INFO.value = -2;
  } else if (M < 0) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (K < 0 || K > Q) {
    INFO.value = -5;
  } else if (NB < 1 || (NB > K && K > 0)) {
    INFO.value = -6;
  } else if (LDV < max(1, Q)) {
    INFO.value = -8;
  } else if (LDT < NB) {
    INFO.value = -10;
  } else if (LDC < max(1, M)) {
    INFO.value = -12;
  }

  if (INFO.value != 0) {
    xerbla('DGEMQRT', -INFO.value);
    return;
  }

  // .. Quick return if possible ..

  if (M == 0 || N == 0 || K == 0) return;

  if (LEFT && TRAN) {
    for (I = 1; NB < 0 ? I >= K : I <= K; I += NB) {
      IB = min(NB, K - I + 1);
      dlarfb('L', 'T', 'F', 'C', M - I + 1, N, IB, V(I, I), LDV, T(1, I), LDT,
          C(I, 1), LDC, WORK.asMatrix(LDWORK), LDWORK);
    }
  } else if (RIGHT && NOTRAN) {
    for (I = 1; NB < 0 ? I >= K : I <= K; I += NB) {
      IB = min(NB, K - I + 1);
      dlarfb('R', 'N', 'F', 'C', M, N - I + 1, IB, V(I, I), LDV, T(1, I), LDT,
          C(1, I), LDC, WORK.asMatrix(LDWORK), LDWORK);
    }
  } else if (LEFT && NOTRAN) {
    KF = ((K - 1) ~/ NB) * NB + 1;
    for (I = KF; -NB < 0 ? I >= 1 : I <= 1; I += -NB) {
      IB = min(NB, K - I + 1);
      dlarfb('L', 'N', 'F', 'C', M - I + 1, N, IB, V(I, I), LDV, T(1, I), LDT,
          C(I, 1), LDC, WORK.asMatrix(LDWORK), LDWORK);
    }
  } else if (RIGHT && TRAN) {
    KF = ((K - 1) ~/ NB) * NB + 1;
    for (I = KF; -NB < 0 ? I >= 1 : I <= 1; I += -NB) {
      IB = min(NB, K - I + 1);
      dlarfb('R', 'T', 'F', 'C', M, N - I + 1, IB, V(I, I), LDV, T(1, I), LDT,
          C(1, I), LDC, WORK.asMatrix(LDWORK), LDWORK);
    }
  }
}
