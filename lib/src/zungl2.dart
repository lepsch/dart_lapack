import 'dart:math';

import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacgv.dart';
import 'package:lapack/src/zlarf.dart';

void zungl2(
  final int M,
  final int N,
  final int K,
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
  int I, J, L;

  // Test the input arguments

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < M) {
    INFO.value = -2;
  } else if (K < 0 || K > M) {
    INFO.value = -3;
  } else if (LDA < max(1, M)) {
    INFO.value = -5;
  }
  if (INFO.value != 0) {
    xerbla('ZUNGL2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (M <= 0) return;

  if (K < M) {
    // Initialise rows k+1:m to rows of the unit matrix

    for (J = 1; J <= N; J++) {
      for (L = K + 1; L <= M; L++) {
        A[L][J] = Complex.zero;
      }
      if (J > K && J <= M) A[J][J] = Complex.one;
    }
  }

  for (I = K; I >= 1; I--) {
    // Apply H(i)**H to A(i:m,i:n) from the right

    if (I < N) {
      zlacgv(N - I, A(I, I + 1).asArray(), LDA);
      if (I < M) {
        A[I][I] = Complex.one;
        zlarf('Right', M - I, N - I + 1, A(I, I).asArray(), LDA,
            TAU[I].conjugate(), A(I + 1, I), LDA, WORK);
      }
      zscal(N - I, -TAU[I], A(I, I + 1).asArray(), LDA);
      zlacgv(N - I, A(I, I + 1).asArray(), LDA);
    }
    A[I][I] = Complex.one - TAU[I].conjugate();

    // Set A(i,1:i-1) to zero

    for (L = 1; L <= I - 1; L++) {
      A[I][L] = Complex.zero;
    }
  }
}
