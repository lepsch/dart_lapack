import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlarf.dart';
import 'package:lapack/src/dlarfg.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgeqr2(
  final int M,
  final int N,
  final Matrix<double> A,
  final int LDA,
  final Array<double> TAU,
  final Array<double> WORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0;
  int I, K;
  double AII;

  // Test the input arguments

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, M)) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('DGEQR2', -INFO.value);
    return;
  }

  K = min(M, N);

  for (I = 1; I <= K; I++) {
    // Generate elementary reflector H(i) to annihilate A(i+1:m,i)

    dlarfg(M - I + 1, A(I, I), A(min(I + 1, M), I), 1, TAU(I));
    if (I < N) {
      // Apply H(i) to A(i:m,i+1:n) from the left

      AII = A[I][I];
      A[I][I] = ONE;
      dlarf(
        'Left',
        M - I + 1,
        N - I,
        A(I, I),
        1,
        TAU(I),
        A(I, I + 1),
        LDA,
        WORK,
      );
      A[I][I] = AII;
    }
  }
}
