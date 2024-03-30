import 'dart:math';

import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dsymv.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/dlarfg.dart';
import 'package:lapack/src/matrix.dart';

void dlatrd(
  final String UPLO,
  final int N,
  final int NB,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> E_,
  final Array<double> TAU_,
  final Matrix<double> W_,
  final int LDW,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final E = E_.having();
  final TAU = TAU_.having();
  final W = W_.having(ld: LDW);
  const ZERO = 0.0, ONE = 1.0, HALF = 0.5;
  int I, IW;
  double ALPHA;

  // Quick return if possible

  if (N <= 0) return;

  if (lsame(UPLO, 'U')) {
    // Reduce last NB columns of upper triangle

    for (I = N; I >= N - NB + 1; I--) {
      IW = I - N + NB;
      if (I < N) {
        // Update A[1:i][i]

        dgemv('No transpose', I, N - I, -ONE, A(1, I + 1), LDA,
            W(I, IW + 1).asArray(), LDW, ONE, A(1, I).asArray(), 1);
        dgemv('No transpose', I, N - I, -ONE, W(1, IW + 1), LDW,
            A(I, I + 1).asArray(), LDA, ONE, A(1, I).asArray(), 1);
      }
      if (I > 1) {
        // Generate elementary reflector H(i) to annihilate
        // A[1:i-2][i]

        dlarfg(I - 1, A.box(I - 1, I), A(1, I).asArray(), 1, TAU.box(I - 1));
        E[I - 1] = A[I - 1][I];
        A[I - 1][I] = ONE;

        // Compute W[1:i-1][i]

        dsymv('Upper', I - 1, ONE, A, LDA, A(1, I).asArray(), 1, ZERO,
            W(1, IW).asArray(), 1);
        if (I < N) {
          dgemv('Transpose', I - 1, N - I, ONE, W(1, IW + 1), LDW,
              A(1, I).asArray(), 1, ZERO, W(I + 1, IW).asArray(), 1);
          dgemv('No transpose', I - 1, N - I, -ONE, A(1, I + 1), LDA,
              W(I + 1, IW).asArray(), 1, ONE, W(1, IW).asArray(), 1);
          dgemv('Transpose', I - 1, N - I, ONE, A(1, I + 1), LDA,
              A(1, I).asArray(), 1, ZERO, W(I + 1, IW).asArray(), 1);
          dgemv('No transpose', I - 1, N - I, -ONE, W(1, IW + 1), LDW,
              W(I + 1, IW).asArray(), 1, ONE, W(1, IW).asArray(), 1);
        }
        dscal(I - 1, TAU[I - 1], W(1, IW).asArray(), 1);
        ALPHA = -HALF *
            TAU[I - 1] *
            ddot(I - 1, W(1, IW).asArray(), 1, A(1, I).asArray(), 1);
        daxpy(I - 1, ALPHA, A(1, I).asArray(), 1, W(1, IW).asArray(), 1);
      }
    }
  } else {
    // Reduce first NB columns of lower triangle

    for (I = 1; I <= NB; I++) {
      // Update A[i:n][i]

      dgemv('No transpose', N - I + 1, I - 1, -ONE, A(I, 1), LDA,
          W(I, 1).asArray(), LDW, ONE, A(I, I).asArray(), 1);
      dgemv('No transpose', N - I + 1, I - 1, -ONE, W(I, 1), LDW,
          A(I, 1).asArray(), LDA, ONE, A(I, I).asArray(), 1);
      if (I < N) {
        // Generate elementary reflector H(i) to annihilate
        // A[i+2:n][i]

        dlarfg(N - I, A.box(I + 1, I), A(min(I + 2, N), I).asArray(), 1,
            TAU.box(I));
        E[I] = A[I + 1][I];
        A[I + 1][I] = ONE;

        // Compute W[i+1:n][i]

        dsymv('Lower', N - I, ONE, A(I + 1, I + 1), LDA, A(I + 1, I).asArray(),
            1, ZERO, W(I + 1, I).asArray(), 1);
        dgemv('Transpose', N - I, I - 1, ONE, W(I + 1, 1), LDW,
            A(I + 1, I).asArray(), 1, ZERO, W(1, I).asArray(), 1);
        dgemv('No transpose', N - I, I - 1, -ONE, A(I + 1, 1), LDA,
            W(1, I).asArray(), 1, ONE, W(I + 1, I).asArray(), 1);
        dgemv('Transpose', N - I, I - 1, ONE, A(I + 1, 1), LDA,
            A(I + 1, I).asArray(), 1, ZERO, W(1, I).asArray(), 1);
        dgemv('No transpose', N - I, I - 1, -ONE, W(I + 1, 1), LDW,
            W(1, I).asArray(), 1, ONE, W(I + 1, I).asArray(), 1);
        dscal(N - I, TAU[I], W(I + 1, I).asArray(), 1);
        ALPHA = -HALF *
            TAU[I] *
            ddot(N - I, W(I + 1, I).asArray(), 1, A(I + 1, I).asArray(), 1);
        daxpy(N - I, ALPHA, A(I + 1, I).asArray(), 1, W(I + 1, I).asArray(), 1);
      }
    }
  }
}
