import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgemlqt.dart';
import 'package:lapack/src/dtpmlqt.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlamswlq(
  final String SIDE,
  final String TRANS,
  final int M,
  final int N,
  final int K,
  final int MB,
  final int NB,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> T_,
  final int LDT,
  final Matrix<double> C_,
  final int LDC,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final T = T_.having(ld: LDT);
  final C = C_.having(ld: LDC);
  final WORK = WORK_.having();
  bool LEFT, RIGHT, TRAN, NOTRAN, LQUERY;
  int I, II, KK, CTR, LW, MINMNK, LWMIN;

  // Test the input arguments

  LQUERY = (LWORK == -1);
  NOTRAN = lsame(TRANS, 'N');
  TRAN = lsame(TRANS, 'T');
  LEFT = lsame(SIDE, 'L');
  RIGHT = lsame(SIDE, 'R');
  if (LEFT) {
    LW = N * MB;
  } else {
    LW = M * MB;
  }

  MINMNK = min(M, min(N, K));
  if (MINMNK == 0) {
    LWMIN = 1;
  } else {
    LWMIN = max(1, LW);
  }

  INFO.value = 0;
  if (!LEFT && !RIGHT) {
    INFO.value = -1;
  } else if (!TRAN && !NOTRAN) {
    INFO.value = -2;
  } else if (K < 0) {
    INFO.value = -5;
  } else if (M < K) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (K < MB || MB < 1) {
    INFO.value = -6;
  } else if (LDA < max(1, K)) {
    INFO.value = -9;
  } else if (LDT < max(1, MB)) {
    INFO.value = -11;
  } else if (LDC < max(1, M)) {
    INFO.value = -13;
  } else if (LWORK < LWMIN && (!LQUERY)) {
    INFO.value = -15;
  }

  if (INFO.value == 0) {
    WORK[1] = LWMIN.toDouble();
  }
  if (INFO.value != 0) {
    xerbla('DLAMSWLQ', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (MINMNK == 0) {
    return;
  }

  if ((NB <= K) || (NB >= max(M, max(N, K)))) {
    dgemlqt(SIDE, TRANS, M, N, K, MB, A, LDA, T, LDT, C, LDC, WORK, INFO);
    return;
  }

  if (LEFT && TRAN) {
    // Multiply Q to the last block of C

    KK = (M - K) % (NB - K);
    CTR = (M - K) ~/ (NB - K);
    if (KK > 0) {
      II = M - KK + 1;
      dtpmlqt('L', 'T', KK, N, K, 0, MB, A(1, II), LDA, T(1, CTR * K + 1), LDT,
          C(1, 1), LDC, C(II, 1), LDC, WORK, INFO);
    } else {
      II = M + 1;
    }

    for (I = II - (NB - K); I >= NB + 1; I -= (NB - K)) {
      // Multiply Q to the current block of C (1:M,I:I+NB)

      CTR = CTR - 1;
      dtpmlqt('L', 'T', NB - K, N, K, 0, MB, A(1, I), LDA, T(1, CTR * K + 1),
          LDT, C(1, 1), LDC, C(I, 1), LDC, WORK, INFO);
    }

    // Multiply Q to the first block of C (1:M,1:NB)

    dgemlqt(
        'L', 'T', NB, N, K, MB, A(1, 1), LDA, T, LDT, C(1, 1), LDC, WORK, INFO);
  } else if (LEFT && NOTRAN) {
    // Multiply Q to the first block of C

    KK = ((M - K) % (NB - K));
    II = M - KK + 1;
    CTR = 1;
    dgemlqt(
        'L', 'N', NB, N, K, MB, A(1, 1), LDA, T, LDT, C(1, 1), LDC, WORK, INFO);

    for (I = NB + 1;
        (NB - K) < 0 ? I >= II - NB + K : I <= II - NB + K;
        I += (NB - K)) {
      // Multiply Q to the current block of C (I:I+NB,1:N)

      dtpmlqt('L', 'N', NB - K, N, K, 0, MB, A(1, I), LDA, T(1, CTR * K + 1),
          LDT, C(1, 1), LDC, C(I, 1), LDC, WORK, INFO);
      CTR = CTR + 1;
    }
    if (II <= M) {
      // Multiply Q to the last block of C

      dtpmlqt('L', 'N', KK, N, K, 0, MB, A(1, II), LDA, T(1, CTR * K + 1), LDT,
          C(1, 1), LDC, C(II, 1), LDC, WORK, INFO);
    }
  } else if (RIGHT && NOTRAN) {
    // Multiply Q to the last block of C

    KK = (N - K) % (NB - K);
    CTR = (N - K) ~/ (NB - K);
    if (KK > 0) {
      II = N - KK + 1;
      dtpmlqt('R', 'N', M, KK, K, 0, MB, A(1, II), LDA, T(1, CTR * K + 1), LDT,
          C(1, 1), LDC, C(1, II), LDC, WORK, INFO);
    } else {
      II = N + 1;
    }

    for (I = II - (NB - K); I >= NB + 1; I -= (NB - K)) {
      // Multiply Q to the current block of C (1:M,I:I+MB)

      CTR = CTR - 1;
      dtpmlqt('R', 'N', M, NB - K, K, 0, MB, A(1, I), LDA, T(1, CTR * K + 1),
          LDT, C(1, 1), LDC, C(1, I), LDC, WORK, INFO);
    }

    // Multiply Q to the first block of C (1:M,1:MB)

    dgemlqt(
        'R', 'N', M, NB, K, MB, A(1, 1), LDA, T, LDT, C(1, 1), LDC, WORK, INFO);
  } else if (RIGHT && TRAN) {
    // Multiply Q to the first block of C

    KK = ((N - K) % (NB - K));
    CTR = 1;
    II = N - KK + 1;
    dgemlqt(
        'R', 'T', M, NB, K, MB, A(1, 1), LDA, T, LDT, C(1, 1), LDC, WORK, INFO);

    for (I = NB + 1;
        (NB - K) < 0 ? I >= II - NB + K : I <= II - NB + K;
        I += (NB - K)) {
      // Multiply Q to the current block of C (1:M,I:I+MB)

      dtpmlqt('R', 'T', M, NB - K, K, 0, MB, A(1, I), LDA, T(1, CTR * K + 1),
          LDT, C(1, 1), LDC, C(1, I), LDC, WORK, INFO);
      CTR = CTR + 1;
    }
    if (II <= N) {
      // Multiply Q to the last block of C

      dtpmlqt('R', 'T', M, KK, K, 0, MB, A(1, II), LDA, T(1, CTR * K + 1), LDT,
          C(1, 1), LDC, C(1, II), LDC, WORK, INFO);
    }
  }

  WORK[1] = LWMIN.toDouble();
}
