import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/ztprfb.dart';

void ztpmlqt(
  final String SIDE,
  final String TRANS,
  final int M,
  final int N,
  final int K,
  final int L,
  final int MB,
  final Matrix<Complex> V_,
  final int LDV,
  final Matrix<Complex> T_,
  final int LDT,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final V = V_.dim(LDV);
  final T = T_.dim(LDT);
  final A = A_.dim(LDA);
  final B = B_.dim(LDB);
  final WORK = WORK_.dim();
  bool LEFT, RIGHT, TRAN, NOTRAN;
  int I, IB, NB, LB, KF, LDAQ = 0;

  // .. Test the input arguments ..

  INFO.value = 0;
  LEFT = lsame(SIDE, 'L');
  RIGHT = lsame(SIDE, 'R');
  TRAN = lsame(TRANS, 'C');
  NOTRAN = lsame(TRANS, 'N');

  if (LEFT) {
    LDAQ = max(1, K);
  } else if (RIGHT) {
    LDAQ = max(1, M);
  }
  if (!LEFT && !RIGHT) {
    INFO.value = -1;
  } else if (!TRAN && !NOTRAN) {
    INFO.value = -2;
  } else if (M < 0) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (K < 0) {
    INFO.value = -5;
  } else if (L < 0 || L > K) {
    INFO.value = -6;
  } else if (MB < 1 || (MB > K && K > 0)) {
    INFO.value = -7;
  } else if (LDV < K) {
    INFO.value = -9;
  } else if (LDT < MB) {
    INFO.value = -11;
  } else if (LDA < LDAQ) {
    INFO.value = -13;
  } else if (LDB < max(1, M)) {
    INFO.value = -15;
  }

  if (INFO.value != 0) {
    xerbla('ZTPMLQT', -INFO.value);
    return;
  }

  // .. Quick return if possible ..

  if (M == 0 || N == 0 || K == 0) return;

  if (LEFT && NOTRAN) {
    for (I = 1; MB < 0 ? I >= K : I <= K; I += MB) {
      IB = min(MB, K - I + 1);
      NB = min(M - L + I + IB - 1, M);
      if (I >= L) {
        LB = 0;
      } else {
        LB = 0;
      }
      ztprfb('L', 'C', 'F', 'R', NB, N, IB, LB, V(I, 1), LDV, T(1, I), LDT,
          A(I, 1), LDA, B, LDB, WORK.asMatrix(IB), IB);
    }
  } else if (RIGHT && TRAN) {
    for (I = 1; MB < 0 ? I >= K : I <= K; I += MB) {
      IB = min(MB, K - I + 1);
      NB = min(N - L + I + IB - 1, N);
      if (I >= L) {
        LB = 0;
      } else {
        LB = NB - N + L - I + 1;
      }
      ztprfb('R', 'N', 'F', 'R', M, NB, IB, LB, V(I, 1), LDV, T(1, I), LDT,
          A(1, I), LDA, B, LDB, WORK.asMatrix(M), M);
    }
  } else if (LEFT && TRAN) {
    KF = ((K - 1) ~/ MB) * MB + 1;
    for (I = KF; -MB < 0 ? I >= 1 : I <= 1; I += -MB) {
      IB = min(MB, K - I + 1);
      NB = min(M - L + I + IB - 1, M);
      if (I >= L) {
        LB = 0;
      } else {
        LB = 0;
      }
      ztprfb('L', 'N', 'F', 'R', NB, N, IB, LB, V(I, 1), LDV, T(1, I), LDT,
          A(I, 1), LDA, B, LDB, WORK.asMatrix(IB), IB);
    }
  } else if (RIGHT && NOTRAN) {
    KF = ((K - 1) ~/ MB) * MB + 1;
    for (I = KF; -MB < 0 ? I >= 1 : I <= 1; I += -MB) {
      IB = min(MB, K - I + 1);
      NB = min(N - L + I + IB - 1, N);
      if (I >= L) {
        LB = 0;
      } else {
        LB = NB - N + L - I + 1;
      }
      ztprfb('R', 'C', 'F', 'R', M, NB, IB, LB, V(I, 1), LDV, T(1, I), LDT,
          A(1, I), LDA, B, LDB, WORK.asMatrix(M), M);
    }
  }
}
