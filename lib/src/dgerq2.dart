import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlarf.dart';
import 'package:lapack/src/dlarfg.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgerq2(
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
  // ..
  // .. External Subroutines ..
  // EXTERNAL DLARF, DLARFG, XERBLA
  // ..
  // .. Intrinsic Functions ..
  // INTRINSIC MAX, MIN

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
    xerbla('DGERQ2', -INFO.value);
    return;
  }

  K = min(M, N);

  for (I = K; I >= 1; I--) {
    // Generate elementary reflector H(i) to annihilate
    // A(m-k+i,1:n-k+i-1)

    dlarfg(N - K + I, A.box(M - K + I, N - K + I), A(M - K + I, 1).asArray(),
        LDA, TAU.box(I));

    // Apply H(i) to A(1:m-k+i-1,1:n-k+i) from the right

    AII = A[M - K + I][N - K + I];
    A[M - K + I][N - K + I] = ONE;
    dlarf('Right', M - K + I - 1, N - K + I, A(M - K + I, 1).asArray(), LDA,
        TAU[I], A, LDA, WORK);
    A[M - K + I][N - K + I] = AII;
  }
}
