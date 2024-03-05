import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgelqt.dart';
import 'package:lapack/src/zlaswlq.dart';

void zgelq(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> T_,
  final int TSIZE,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --
  final A = A_.having(ld: LDA);
  final T = T_.having();
  final WORK = WORK_.having();
  bool LQUERY, LMINWS, MINT, MINW;
  int MB, NB, MINTSZ, NBLCKS, LWMIN, LWOPT, LWREQ;

  // Test the input arguments

  INFO.value = 0;

  LQUERY = (TSIZE == -1 || TSIZE == -2 || LWORK == -1 || LWORK == -2);

  MINT = false;
  MINW = false;
  if (TSIZE == -2 || LWORK == -2) {
    if (TSIZE != -1) MINT = true;
    if (LWORK != -1) MINW = true;
  }

  // Determine the block size

  if (min(M, N) > 0) {
    MB = ilaenv(1, 'ZGELQ ', ' ', M, N, 1, -1);
    NB = ilaenv(1, 'ZGELQ ', ' ', M, N, 2, -1);
  } else {
    MB = 1;
    NB = N;
  }
  if (MB > min(M, N) || MB < 1) MB = 1;
  if (NB > N || NB <= M) NB = N;
  MINTSZ = M + 5;
  if (NB > M && N > M) {
    if (((N - M) % (NB - M)) == 0) {
      NBLCKS = (N - M) ~/ (NB - M);
    } else {
      NBLCKS = (N - M) ~/ (NB - M) + 1;
    }
  } else {
    NBLCKS = 1;
  }

  // Determine if the workspace size satisfies minimal size

  if ((N <= M) || (NB <= M) || (NB >= N)) {
    LWMIN = max(1, N);
    LWOPT = max(1, MB * N);
  } else {
    LWMIN = max(1, M);
    LWOPT = max(1, MB * M);
  }
  LMINWS = false;
  if ((TSIZE < max(1, MB * M * NBLCKS + 5) || LWORK < LWOPT) &&
      (LWORK >= LWMIN) &&
      (TSIZE >= MINTSZ) &&
      (!LQUERY)) {
    if (TSIZE < max(1, MB * M * NBLCKS + 5)) {
      LMINWS = true;
      MB = 1;
      NB = N;
    }
    if (LWORK < LWOPT) {
      LMINWS = true;
      MB = 1;
    }
  }
  if ((N <= M) || (NB <= M) || (NB >= N)) {
    LWREQ = max(1, MB * N);
  } else {
    LWREQ = max(1, MB * M);
  }

  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, M)) {
    INFO.value = -4;
  } else if (TSIZE < max(1, MB * M * NBLCKS + 5) && (!LQUERY) && (!LMINWS)) {
    INFO.value = -6;
  } else if ((LWORK < LWREQ) && (!LQUERY) && (!LMINWS)) {
    INFO.value = -8;
  }

  if (INFO.value == 0) {
    if (MINT) {
      T[1] = MINTSZ.toComplex();
    } else {
      T[1] = (MB * M * NBLCKS + 5).toComplex();
    }
    T[2] = MB.toComplex();
    T[3] = NB.toComplex();
    if (MINW) {
      WORK[1] = LWMIN.toComplex();
    } else {
      WORK[1] = LWREQ.toComplex();
    }
  }
  if (INFO.value != 0) {
    xerbla('ZGELQ', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (min(M, N) == 0) {
    return;
  }

  // The LQ Decomposition

  if ((N <= M) || (NB <= M) || (NB >= N)) {
    zgelqt(M, N, MB, A, LDA, T(6).asMatrix(MB), MB, WORK, INFO);
  } else {
    zlaswlq(M, N, MB, NB, A, LDA, T(6).asMatrix(MB), MB, WORK, LWORK, INFO);
  }

  WORK[1] = LWREQ.toComplex();
}
