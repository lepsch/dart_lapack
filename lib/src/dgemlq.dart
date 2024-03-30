import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgemlqt.dart';
import 'package:lapack/src/dlamswlq.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgemlq(
  final String SIDE,
  final String TRANS,
  final int M,
  final int N,
  final int K,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> T_,
  final int TSIZE,
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
  final T = T_.having();
  final C = C_.having(ld: LDC);
  final WORK = WORK_.having();
  bool LEFT, RIGHT, TRAN, NOTRAN, LQUERY;
  int MB,
      NB,
      LW,
      // NBLCKS,
      MN,
      MINMNK,
      LWMIN;

  // Test the input arguments

  LQUERY = (LWORK == -1);
  NOTRAN = lsame(TRANS, 'N');
  TRAN = lsame(TRANS, 'T');
  LEFT = lsame(SIDE, 'L');
  RIGHT = lsame(SIDE, 'R');

  MB = T[2].toInt();
  NB = T[3].toInt();
  if (LEFT) {
    LW = N * MB;
    MN = M;
  } else {
    LW = M * MB;
    MN = N;
  }

  MINMNK = min(M, min(N, K));
  if (MINMNK == 0) {
    LWMIN = 1;
  } else {
    LWMIN = max(1, LW);
  }

  // if ((NB > K) && (MN > K)) {
  //   if (((MN - K) % (NB - K)) == 0) {
  //     NBLCKS = (MN - K) ~/ (NB - K);
  //   } else {
  //     NBLCKS = (MN - K) ~/ (NB - K) + 1;
  //   }
  // } else {
  //   NBLCKS = 1;
  // }

  INFO.value = 0;
  if (!LEFT && !RIGHT) {
    INFO.value = -1;
  } else if (!TRAN && !NOTRAN) {
    INFO.value = -2;
  } else if (M < 0) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (K < 0 || K > MN) {
    INFO.value = -5;
  } else if (LDA < max(1, K)) {
    INFO.value = -7;
  } else if (TSIZE < 5) {
    INFO.value = -9;
  } else if (LDC < max(1, M)) {
    INFO.value = -11;
  } else if (LWORK < LWMIN && !LQUERY) {
    INFO.value = -13;
  }

  if (INFO.value == 0) {
    WORK[1] = LWMIN.toDouble();
  }

  if (INFO.value != 0) {
    xerbla('DGEMLQ', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (MINMNK == 0) {
    return;
  }

  if ((LEFT && M <= K) ||
      (RIGHT && N <= K) ||
      (NB <= K) ||
      (NB >= max(M, max(N, K)))) {
    dgemlqt(SIDE, TRANS, M, N, K, MB, A, LDA, T(6).asMatrix(MB), MB, C, LDC,
        WORK, INFO);
  } else {
    dlamswlq(SIDE, TRANS, M, N, K, MB, NB, A, LDA, T(6).asMatrix(MB), MB, C,
        LDC, WORK, LWORK, INFO);
  }

  WORK[1] = LWMIN.toDouble();
}
