import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgemlqt.dart';
import 'package:lapack/src/zlamswlq.dart';

void zgemlq(
  final String SIDE,
  final String TRANS,
  final int M,
  final int N,
  final int K,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> T_,
  final int TSIZE,
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
  final C = C_.having(ld: LDC);
  final T = T_.having();
  final WORK = WORK_.having();
  bool LEFT, RIGHT, TRAN, NOTRAN, LQUERY;
  int MB, NB, LW, MN, MINMNK, LWMIN;

  // Test the input arguments

  LQUERY = (LWORK == -1);
  NOTRAN = lsame(TRANS, 'N');
  TRAN = lsame(TRANS, 'C');
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

  // if ( ( NB > K ) && ( MN > K ) ) {
  //   if ( ((MN - K) % (NB - K)) == 0 ) {
  //     NBLCKS = ( MN - K ) ~/ ( NB - K );
  //   } else {
  //     NBLCKS = ( MN - K ) ~/ ( NB - K ) + 1;
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
  } else if ((LWORK < LWMIN) && !LQUERY) {
    INFO.value = -13;
  }

  if (INFO.value == 0) {
    WORK[1] = LW.toComplex();
  }

  if (INFO.value != 0) {
    xerbla('ZGEMLQ', -INFO.value);
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
    zgemlqt(SIDE, TRANS, M, N, K, MB, A, LDA, T(6).asMatrix(MB), MB, C, LDC,
        WORK, INFO);
  } else {
    zlamswlq(SIDE, TRANS, M, N, K, MB, NB, A, LDA, T(6).asMatrix(MB), MB, C,
        LDC, WORK, LWORK, INFO);
  }

  WORK[1] = LW.toComplex();
}
