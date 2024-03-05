import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlarf.dart';
import 'package:lapack/src/dlarfg.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgebd2(
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> D_,
  final Array<double> E_,
  final Array<double> TAUQ_,
  final Array<double> TAUP_,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final D = D_.having();
  final E = E_.having();
  final TAUQ = TAUQ_.having();
  final TAUP = TAUP_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  int I;

  // Test the input parameters

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, M)) {
    INFO.value = -4;
  }
  if (INFO.value < 0) {
    xerbla('DGEBD2', -INFO.value);
    return;
  }

  if (M >= N) {
    // Reduce to upper bidiagonal form

    for (I = 1; I <= N; I++) {
      // Generate elementary reflector H(i) to annihilate A(i+1:m,i)

      dlarfg(M - I + 1, A.box(I, I), A(min(I + 1, M), I).asArray(), 1,
          TAUQ.box(I));
      D[I] = A[I][I];
      A[I][I] = ONE;

      // Apply H(i) to A(i:m,i+1:n) from the left

      if (I < N) {
        dlarf('Left', M - I + 1, N - I, A(I, I).asArray(), 1, TAUQ[I],
            A(I, I + 1), LDA, WORK);
      }
      A[I][I] = D[I];

      if (I < N) {
        // Generate elementary reflector G(i) to annihilate
        // A(i,i+2:n)

        dlarfg(N - I, A.box(I, I + 1), A(I, min(I + 2, N)).asArray(), LDA,
            TAUP.box(I));
        E[I] = A[I][I + 1];
        A[I][I + 1] = ONE;

        // Apply G(i) to A(i+1:m,i+1:n) from the right

        dlarf('Right', M - I, N - I, A(I, I + 1).asArray(), LDA, TAUP[I],
            A(I + 1, I + 1), LDA, WORK);
        A[I][I + 1] = E[I];
      } else {
        TAUP[I] = ZERO;
      }
    }
  } else {
    // Reduce to lower bidiagonal form

    for (I = 1; I <= M; I++) {
      // Generate elementary reflector G(i) to annihilate A(i,i+1:n)

      dlarfg(N - I + 1, A.box(I, I), A(I, min(I + 1, N)).asArray(), LDA,
          TAUP.box(I));
      D[I] = A[I][I];
      A[I][I] = ONE;

      // Apply G(i) to A(i+1:m,i:n) from the right

      if (I < M) {
        dlarf('Right', M - I, N - I + 1, A(I, I).asArray(), LDA, TAUP[I],
            A(I + 1, I), LDA, WORK);
      }
      A[I][I] = D[I];

      if (I < M) {
        // Generate elementary reflector H(i) to annihilate
        // A(i+2:m,i)

        dlarfg(M - I, A.box(I + 1, I), A(min(I + 2, M), I).asArray(), 1,
            TAUQ.box(I));
        E[I] = A[I + 1][I];
        A[I + 1][I] = ONE;

        // Apply H(i) to A(i+1:m,i+1:n) from the left

        dlarf('Left', M - I, N - I, A(I + 1, I).asArray(), 1, TAUQ[I],
            A(I + 1, I + 1), LDA, WORK);
        A[I + 1][I] = E[I];
      } else {
        TAUQ[I] = ZERO;
      }
    }
  }
}
