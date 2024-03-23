import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dtprfb.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dtpmqrt(
  final String SIDE,
  final String TRANS,
  final int M,
  final int N,
  final int K,
  final int L,
  final int NB,
  final Matrix<double> V_,
  final int LDV,
  final Matrix<double> T_,
  final int LDT,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final V = V_.having(ld: LDV);
  final T = T_.having(ld: LDT);
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final WORK = WORK_.having();
  bool LEFT, RIGHT, TRAN, NOTRAN;
  int I, IB, MB, LB, KF, LDAQ = 0, LDVQ = 0;

  // .. Test the input arguments ..

  INFO.value = 0;
  LEFT = lsame(SIDE, 'L');
  RIGHT = lsame(SIDE, 'R');
  TRAN = lsame(TRANS, 'T');
  NOTRAN = lsame(TRANS, 'N');

  if (LEFT) {
    LDVQ = max(1, M);
    LDAQ = max(1, K);
  } else if (RIGHT) {
    LDVQ = max(1, N);
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
  } else if (NB < 1 || (NB > K && K > 0)) {
    INFO.value = -7;
  } else if (LDV < LDVQ) {
    INFO.value = -9;
  } else if (LDT < NB) {
    INFO.value = -11;
  } else if (LDA < LDAQ) {
    INFO.value = -13;
  } else if (LDB < max(1, M)) {
    INFO.value = -15;
  }

  if (INFO.value != 0) {
    xerbla('DTPMQRT', -INFO.value);
    return;
  }

  // .. Quick return if possible ..

  if (M == 0 || N == 0 || K == 0) return;

  if (LEFT && TRAN) {
    for (I = 1; I <= K; I += NB) {
      IB = min(NB, K - I + 1);
      MB = min(M - L + I + IB - 1, M);
      if (I >= L) {
        LB = 0;
      } else {
        LB = MB - M + L - I + 1;
      }
      dtprfb('L', 'T', 'F', 'C', MB, N, IB, LB, V(1, I), LDV, T(1, I), LDT,
          A(I, 1), LDA, B, LDB, WORK.asMatrix(IB), IB);
    }
  } else if (RIGHT && NOTRAN) {
    for (I = 1; I <= K; I += NB) {
      IB = min(NB, K - I + 1);
      MB = min(N - L + I + IB - 1, N);
      if (I >= L) {
        LB = 0;
      } else {
        LB = MB - N + L - I + 1;
      }
      dtprfb('R', 'N', 'F', 'C', M, MB, IB, LB, V(1, I), LDV, T(1, I), LDT,
          A(1, I), LDA, B, LDB, WORK.asMatrix(M), M);
    }
  } else if (LEFT && NOTRAN) {
    KF = ((K - 1) ~/ NB) * NB + 1;
    for (I = KF; I >= 1; I -= NB) {
      IB = min(NB, K - I + 1);
      MB = min(M - L + I + IB - 1, M);
      if (I >= L) {
        LB = 0;
      } else {
        LB = MB - M + L - I + 1;
      }
      dtprfb('L', 'N', 'F', 'C', MB, N, IB, LB, V(1, I), LDV, T(1, I), LDT,
          A(I, 1), LDA, B, LDB, WORK.asMatrix(IB), IB);
    }
  } else if (RIGHT && TRAN) {
    KF = ((K - 1) ~/ NB) * NB + 1;
    for (I = KF; I >= 1; I -= NB) {
      IB = min(NB, K - I + 1);
      MB = min(N - L + I + IB - 1, N);
      if (I >= L) {
        LB = 0;
      } else {
        LB = MB - N + L - I + 1;
      }
      dtprfb('R', 'T', 'F', 'C', M, MB, IB, LB, V(1, I), LDV, T(1, I), LDT,
          A(1, I), LDA, B, LDB, WORK.asMatrix(M), M);
    }
  }
}
