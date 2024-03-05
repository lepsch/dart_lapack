import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlarfb.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgemlqt(
  final String SIDE,
  final String TRANS,
  final int M,
  final int N,
  final int K,
  final int MB,
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
  final V = V_.having(ld: LDV);
  final T = T_.having(ld: LDT);
  final C = C_.having(ld: LDC);
  final WORK = WORK_.having();
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
  } else if (MB < 1 || (MB > K && K > 0)) {
    INFO.value = -6;
  } else if (LDV < max(1, K)) {
    INFO.value = -8;
  } else if (LDT < MB) {
    INFO.value = -10;
  } else if (LDC < max(1, M)) {
    INFO.value = -12;
  }

  if (INFO.value != 0) {
    xerbla('DGEMLQT', -INFO.value);
    return;
  }

  // .. Quick return if possible ..

  if (M == 0 || N == 0 || K == 0) return;

  if (LEFT && NOTRAN) {
    for (I = 1; MB < 0 ? I >= K : I <= K; I += MB) {
      IB = min(MB, K - I + 1);
      dlarfb('L', 'T', 'F', 'R', M - I + 1, N, IB, V(I, I), LDV, T(1, I), LDT,
          C(I, 1), LDC, WORK.asMatrix(LDWORK), LDWORK);
    }
  } else if (RIGHT && TRAN) {
    for (I = 1; MB < 0 ? I >= K : I <= K; I += MB) {
      IB = min(MB, K - I + 1);
      dlarfb('R', 'N', 'F', 'R', M, N - I + 1, IB, V(I, I), LDV, T(1, I), LDT,
          C(1, I), LDC, WORK.asMatrix(LDWORK), LDWORK);
    }
  } else if (LEFT && TRAN) {
    KF = ((K - 1) ~/ MB) * MB + 1;
    for (I = KF; -MB < 0 ? I >= 1 : I <= 1; I += -MB) {
      IB = min(MB, K - I + 1);
      dlarfb('L', 'N', 'F', 'R', M - I + 1, N, IB, V(I, I), LDV, T(1, I), LDT,
          C(I, 1), LDC, WORK.asMatrix(LDWORK), LDWORK);
    }
  } else if (RIGHT && NOTRAN) {
    KF = ((K - 1) ~/ MB) * MB + 1;
    for (I = KF; -MB < 0 ? I >= 1 : I <= 1; I += -MB) {
      IB = min(MB, K - I + 1);
      dlarfb('R', 'T', 'F', 'R', M, N - I + 1, IB, V(I, I), LDV, T(1, I), LDT,
          C(1, I), LDC, WORK.asMatrix(LDWORK), LDWORK);
    }
  }
}
