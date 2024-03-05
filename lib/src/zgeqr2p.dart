import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlarf.dart';
import 'package:lapack/src/zlarfgp.dart';

void zgeqr2p(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  int I, K;
  Complex ALPHA;

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
    xerbla('ZGEQR2P', -INFO.value);
    return;
  }

  K = min(M, N);

  for (I = 1; I <= K; I++) {
    // 10

    // Generate elementary reflector H(i) to annihilate A(i+1:m,i)

    zlarfgp(M - I + 1, A(I, I), A(min(I + 1, M), I).asArray(), 1, TAU(I));
    if (I < N) {
      // Apply H(i)**H to A(i:m,i+1:n) from the left

      ALPHA = A[I][I];
      A[I][I] = Complex.one;
      zlarf('Left', M - I + 1, N - I, A(I, I).asArray(), 1, TAU[I].conjugate(),
          A(I, I + 1), LDA, WORK);
      A[I][I] = ALPHA;
    }
  } // 10
}
