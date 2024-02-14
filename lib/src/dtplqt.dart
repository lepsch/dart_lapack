import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dtplqt2.dart';
import 'package:lapack/src/dtprfb.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dtplqt(
  final int M,
  final int N,
  final int L,
  final int MB,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> T_,
  final int LDT,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final B = B_.dim(LDB);
  final T = T_.dim(LDT);
  final WORK = WORK_.dim();
  int I, IB, LB, NB;
  final IINFO = Box(0);

  // Test the input arguments

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (L < 0 || (L > min(M, N) && min(M, N) >= 0)) {
    INFO.value = -3;
  } else if (MB < 1 || (MB > M && M > 0)) {
    INFO.value = -4;
  } else if (LDA < max(1, M)) {
    INFO.value = -6;
  } else if (LDB < max(1, M)) {
    INFO.value = -8;
  } else if (LDT < MB) {
    INFO.value = -10;
  }
  if (INFO.value != 0) {
    xerbla('DTPLQT', -INFO.value);
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) return;

  for (I = 1; MB < 0 ? I >= M : I <= M; I += MB) {
    // Compute the QR factorization of the current block

    IB = min(M - I + 1, MB);
    NB = min(N - L + I + IB - 1, N);
    if (I >= L) {
      LB = 0;
    } else {
      LB = NB - N + L - I + 1;
    }

    dtplqt2(IB, NB, LB, A(I, I), LDA, B(I, 1), LDB, T(1, I), LDT, IINFO);

    // Update by applying H**T to B(I+IB:M,:) from the right

    if (I + IB <= M) {
      dtprfb(
        'R',
        'N',
        'F',
        'R',
        M - I - IB + 1,
        NB,
        IB,
        LB,
        B(I, 1),
        LDB,
        T(1, I),
        LDT,
        A(I + IB, I),
        LDA,
        B(I + IB, 1),
        LDB,
        WORK.asMatrix(M - I - IB + 1),
        M - I - IB + 1,
      );
    }
  }
}
