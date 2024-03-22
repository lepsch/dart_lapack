import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgemqrt.dart';
import 'package:lapack/src/ztpmqrt.dart';

void zlamtsqr(
  final String SIDE,
  final String TRANS,
  final int M,
  final int N,
  final int K,
  final int MB,
  final int NB,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> T_,
  final int LDT,
  final Matrix<Complex> C_,
  final int LDC,
  final Array<Complex> WORK_,
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
  int I, II, KK, LW, CTR, Q, MINMNK, LWMIN;

  // Test the input arguments

  INFO.value = 0;
  LQUERY = (LWORK == -1);
  NOTRAN = lsame(TRANS, 'N');
  TRAN = lsame(TRANS, 'C');
  LEFT = lsame(SIDE, 'L');
  RIGHT = lsame(SIDE, 'R');
  if (LEFT) {
    LW = N * NB;
    Q = M;
  } else {
    LW = M * NB;
    Q = N;
  }

  MINMNK = min(M, min(N, K));
  if (MINMNK == 0) {
    LWMIN = 1;
  } else {
    LWMIN = max(1, LW);
  }

  if (!LEFT && !RIGHT) {
    INFO.value = -1;
  } else if (!TRAN && !NOTRAN) {
    INFO.value = -2;
  } else if (M < K) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (K < 0) {
    INFO.value = -5;
  } else if (K < NB || NB < 1) {
    INFO.value = -7;
  } else if (LDA < max(1, Q)) {
    INFO.value = -9;
  } else if (LDT < max(1, NB)) {
    INFO.value = -11;
  } else if (LDC < max(1, M)) {
    INFO.value = -13;
  } else if (LWORK < LWMIN && !LQUERY) {
    INFO.value = -15;
  }

  if (INFO.value == 0) {
    WORK[1] = LWMIN.toComplex();
  }

  if (INFO.value != 0) {
    xerbla('ZLAMTSQR', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (MINMNK == 0) {
    return;
  }

  // Determine the block size if it is tall skinny or short and wide

  if ((MB <= K) || (MB >= max(M, min(N, K)))) {
    zgemqrt(SIDE, TRANS, M, N, K, NB, A, LDA, T, LDT, C, LDC, WORK, INFO);
    return;
  }

  if (LEFT && NOTRAN) {
    // Multiply Q to the last block of C

    KK = ((M - K) % (MB - K));
    CTR = (M - K) ~/ (MB - K);
    if (KK > 0) {
      II = M - KK + 1;
      ztpmqrt('L', 'N', KK, N, K, 0, NB, A(II, 1), LDA, T(1, CTR * K + 1), LDT,
          C(1, 1), LDC, C(II, 1), LDC, WORK, INFO);
    } else {
      II = M + 1;
    }

    for (I = II - (MB - K);
        -(MB - K) < 0 ? I >= MB + 1 : I <= MB + 1;
        I -= (MB - K)) {
      // Multiply Q to the current block of C (I:I+MB,1:N)

      CTR--;
      ztpmqrt('L', 'N', MB - K, N, K, 0, NB, A(I, 1), LDA, T(1, CTR * K + 1),
          LDT, C(1, 1), LDC, C(I, 1), LDC, WORK, INFO);
    }

    // Multiply Q to the first block of C (1:MB,1:N)

    zgemqrt(
        'L', 'N', MB, N, K, NB, A(1, 1), LDA, T, LDT, C(1, 1), LDC, WORK, INFO);
  } else if (LEFT && TRAN) {
    // Multiply Q to the first block of C

    KK = ((M - K) % (MB - K));
    II = M - KK + 1;
    CTR = 1;
    zgemqrt(
        'L', 'C', MB, N, K, NB, A(1, 1), LDA, T, LDT, C(1, 1), LDC, WORK, INFO);

    for (I = MB + 1;
        (MB - K) < 0 ? I >= II - MB + K : I <= II - MB + K;
        I += (MB - K)) {
      // Multiply Q to the current block of C (I:I+MB,1:N)

      ztpmqrt('L', 'C', MB - K, N, K, 0, NB, A(I, 1), LDA, T(1, CTR * K + 1),
          LDT, C(1, 1), LDC, C(I, 1), LDC, WORK, INFO);
      CTR++;
    }
    if (II <= M) {
      // Multiply Q to the last block of C

      ztpmqrt('L', 'C', KK, N, K, 0, NB, A(II, 1), LDA, T(1, CTR * K + 1), LDT,
          C(1, 1), LDC, C(II, 1), LDC, WORK, INFO);
    }
  } else if (RIGHT && TRAN) {
    // Multiply Q to the last block of C

    KK = ((N - K) % (MB - K));
    CTR = (N - K) ~/ (MB - K);
    if (KK > 0) {
      II = N - KK + 1;
      ztpmqrt('R', 'C', M, KK, K, 0, NB, A(II, 1), LDA, T(1, CTR * K + 1), LDT,
          C(1, 1), LDC, C(1, II), LDC, WORK, INFO);
    } else {
      II = N + 1;
    }

    for (I = II - (MB - K);
        -(MB - K) < 0 ? I >= MB + 1 : I <= MB + 1;
        I -= (MB - K)) {
      // Multiply Q to the current block of C (1:M,I:I+MB)

      CTR--;
      ztpmqrt('R', 'C', M, MB - K, K, 0, NB, A(I, 1), LDA, T(1, CTR * K + 1),
          LDT, C(1, 1), LDC, C(1, I), LDC, WORK, INFO);
    }

    // Multiply Q to the first block of C (1:M,1:MB)

    zgemqrt(
        'R', 'C', M, MB, K, NB, A(1, 1), LDA, T, LDT, C(1, 1), LDC, WORK, INFO);
  } else if (RIGHT && NOTRAN) {
    // Multiply Q to the first block of C

    KK = ((N - K) % (MB - K));
    II = N - KK + 1;
    CTR = 1;
    zgemqrt(
        'R', 'N', M, MB, K, NB, A(1, 1), LDA, T, LDT, C(1, 1), LDC, WORK, INFO);

    for (I = MB + 1;
        (MB - K) < 0 ? I >= II - MB + K : I <= II - MB + K;
        I += (MB - K)) {
      // Multiply Q to the current block of C (1:M,I:I+MB)

      ztpmqrt('R', 'N', M, MB - K, K, 0, NB, A(I, 1), LDA, T(1, CTR * K + 1),
          LDT, C(1, 1), LDC, C(1, I), LDC, WORK, INFO);
      CTR++;
    }
    if (II <= N) {
      // Multiply Q to the last block of C

      ztpmqrt('R', 'N', M, KK, K, 0, NB, A(II, 1), LDA, T(1, CTR * K + 1), LDT,
          C(1, 1), LDC, C(1, II), LDC, WORK, INFO);
    }
  }

  WORK[1] = LWMIN.toComplex();
}
