import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgelqt.dart';
import 'package:lapack/src/ztplqt.dart';

void zlaswlq(
  final int M,
  final int N,
  final int MB,
  final int NB,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> T_,
  final int LDT,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --
  final A = A_.having(ld: LDA);
  final T = T_.having(ld: LDT);
  final WORK = WORK_.having();
  bool LQUERY;
  int I, II, KK, CTR, MINMN, LWMIN;

  // TEST THE INPUT ARGUMENTS

  INFO.value = 0;

  LQUERY = (LWORK == -1);

  MINMN = min(M, N);
  if (MINMN == 0) {
    LWMIN = 1;
  } else {
    LWMIN = M * MB;
  }

  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0 || N < M) {
    INFO.value = -2;
  } else if (MB < 1 || (MB > M && M > 0)) {
    INFO.value = -3;
  } else if (NB <= 0) {
    INFO.value = -4;
  } else if (LDA < max(1, M)) {
    INFO.value = -6;
  } else if (LDT < MB) {
    INFO.value = -8;
  } else if (LWORK < LWMIN && !LQUERY) {
    INFO.value = -10;
  }

  if (INFO.value == 0) {
    WORK[1] = LWMIN.toComplex();
  }

  if (INFO.value != 0) {
    xerbla('ZLASWLQ', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (MINMN == 0) {
    return;
  }

  // The LQ Decomposition

  if ((M >= N) || (NB <= M) || (NB >= N)) {
    zgelqt(M, N, MB, A, LDA, T, LDT, WORK, INFO);
    return;
  }

  KK = ((N - M) % (NB - M));
  II = N - KK + 1;

  // Compute the LQ factorization of the first block A(1:M,1:NB)

  zgelqt(M, NB, MB, A(1, 1), LDA, T, LDT, WORK, INFO);
  CTR = 1;

  for (I = NB + 1;
      (NB - M) < 0 ? I >= II - NB + M : I <= II - NB + M;
      I += (NB - M)) {
    // Compute the QR factorization of the current block A(1:M,I:I+NB-M)

    ztplqt(M, NB - M, 0, MB, A(1, 1), LDA, A(1, I), LDA, T(1, CTR * M + 1), LDT,
        WORK, INFO);
    CTR++;
  }

  // Compute the QR factorization of the last block A(1:M,II:N)

  if (II <= N) {
    ztplqt(M, KK, 0, MB, A(1, 1), LDA, A(1, II), LDA, T(1, CTR * M + 1), LDT,
        WORK, INFO);
  }

  WORK[1] = LWMIN.toComplex();
}
