import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgelqt3.dart';
import 'package:lapack/src/dlarfb.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgelqt(
  final int M,
  final int N,
  final int MB,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> T_,
  final int LDT,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final T = T_.having(ld: LDT);
  final WORK = WORK_.having();
  int I, IB, K;
  final IINFO = Box(0);

  // Test the input arguments

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (MB < 1 || (MB > min(M, N) && min(M, N) > 0)) {
    INFO.value = -3;
  } else if (LDA < max(1, M)) {
    INFO.value = -5;
  } else if (LDT < MB) {
    INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('DGELQT', -INFO.value);
    return;
  }

  // Quick return if possible

  K = min(M, N);
  if (K == 0) return;

  // Blocked loop of length K

  for (I = 1; MB < 0 ? I >= K : I <= K; I += MB) {
    IB = min(K - I + 1, MB);

    // Compute the LQ factorization of the current block A(I:M,I:I+IB-1)

    dgelqt3(IB, N - I + 1, A(I, I), LDA, T(1, I), LDT, IINFO);
    if (I + IB <= M) {
      // Update by applying H**T to A(I:M,I+IB:N) from the right

      dlarfb(
          'R',
          'N',
          'F',
          'R',
          M - I - IB + 1,
          N - I + 1,
          IB,
          A(I, I),
          LDA,
          T(1, I),
          LDT,
          A(I + IB, I),
          LDA,
          WORK.asMatrix(M - I - IB + 1),
          M - I - IB + 1);
    }
  }
}
