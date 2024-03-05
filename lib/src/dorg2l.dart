import 'dart:math';

import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlarf.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dorg2l(
  final int M,
  final int N,
  final int K,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  int I, II, J, L;

  // Test the input arguments

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0 || N > M) {
    INFO.value = -2;
  } else if (K < 0 || K > N) {
    INFO.value = -3;
  } else if (LDA < max(1, M)) {
    INFO.value = -5;
  }
  if (INFO.value != 0) {
    xerbla('DORG2L', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N <= 0) return;

  // Initialise columns 1:n-k to columns of the unit matrix

  for (J = 1; J <= N - K; J++) {
    for (L = 1; L <= M; L++) {
      A[L][J] = ZERO;
    }
    A[M - N + J][J] = ONE;
  }

  for (I = 1; I <= K; I++) {
    II = N - K + I;

    // Apply H(i) to A(1:m-k+i,1:n-k+i) from the left

    A[M - N + II][II] = ONE;
    dlarf('Left', M - N + II, II - 1, A(1, II).asArray(), 1, TAU[I], A, LDA,
        WORK);
    dscal(M - N + II - 1, -TAU[I], A(1, II).asArray(), 1);
    A[M - N + II][II] = ONE - TAU[I];

    // Set A(m-k+i+1:m,n-k+i) to zero

    for (L = M - N + II + 1; L <= M; L++) {
      A[L][II] = ZERO;
    }
  }
}
