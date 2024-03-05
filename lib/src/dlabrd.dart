import 'dart:math';

import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/dlarfg.dart';
import 'package:lapack/src/matrix.dart';

void dlabrd(
  final int M,
  final int N,
  final int NB,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> D_,
  final Array<double> E_,
  final Array<double> TAUQ_,
  final Array<double> TAUP_,
  final Matrix<double> X_,
  final int LDX,
  final Matrix<double> Y_,
  final int LDY,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final D = D_.having();
  final E = E_.having();
  final TAUQ = TAUQ_.having();
  final TAUP = TAUP_.having();
  final X = X_.having(ld: LDX);
  final Y = Y_.having(ld: LDY);
  const ZERO = 0.0, ONE = 1.0;
  int I;

  // Quick return if possible

  if (M <= 0 || N <= 0) return;

  if (M >= N) {
    // Reduce to upper bidiagonal form

    for (I = 1; I <= NB; I++) {
      // Update A(i:m,i)

      dgemv('No transpose', M - I + 1, I - 1, -ONE, A(I, 1), LDA,
          Y(I, 1).asArray(), LDY, ONE, A(I, I).asArray(), 1);
      dgemv('No transpose', M - I + 1, I - 1, -ONE, X(I, 1), LDX,
          A(1, I).asArray(), 1, ONE, A(I, I).asArray(), 1);

      // Generate reflection Q(i) to annihilate A(i+1:m,i)

      dlarfg(M - I + 1, A.box(I, I), A(min(I + 1, M), I).asArray(), 1,
          TAUQ.box(I));
      D[I] = A[I][I];
      if (I < N) {
        A[I][I] = ONE;

        // Compute Y(i+1:n,i)

        dgemv('Transpose', M - I + 1, N - I, ONE, A(I, I + 1), LDA,
            A(I, I).asArray(), 1, ZERO, Y(I + 1, I).asArray(), 1);
        dgemv('Transpose', M - I + 1, I - 1, ONE, A(I, 1), LDA,
            A(I, I).asArray(), 1, ZERO, Y(1, I).asArray(), 1);
        dgemv('No transpose', N - I, I - 1, -ONE, Y(I + 1, 1), LDY,
            Y(1, I).asArray(), 1, ONE, Y(I + 1, I).asArray(), 1);
        dgemv('Transpose', M - I + 1, I - 1, ONE, X(I, 1), LDX,
            A(I, I).asArray(), 1, ZERO, Y(1, I).asArray(), 1);
        dgemv('Transpose', I - 1, N - I, -ONE, A(1, I + 1), LDA,
            Y(1, I).asArray(), 1, ONE, Y(I + 1, I).asArray(), 1);
        dscal(N - I, TAUQ[I], Y(I + 1, I).asArray(), 1);

        // Update A(i,i+1:n)

        dgemv('No transpose', N - I, I, -ONE, Y(I + 1, 1), LDY,
            A(I, 1).asArray(), LDA, ONE, A(I, I + 1).asArray(), LDA);
        dgemv('Transpose', I - 1, N - I, -ONE, A(1, I + 1), LDA,
            X(I, 1).asArray(), LDX, ONE, A(I, I + 1).asArray(), LDA);

        // Generate reflection P(i) to annihilate A(i,i+2:n)

        dlarfg(N - I, A.box(I, I + 1), A(I, min(I + 2, N)).asArray(), LDA,
            TAUP.box(I));
        E[I] = A[I][I + 1];
        A[I][I + 1] = ONE;

        // Compute X(i+1:m,i)

        dgemv('No transpose', M - I, N - I, ONE, A(I + 1, I + 1), LDA,
            A(I, I + 1).asArray(), LDA, ZERO, X(I + 1, I).asArray(), 1);
        dgemv('Transpose', N - I, I, ONE, Y(I + 1, 1), LDY,
            A(I, I + 1).asArray(), LDA, ZERO, X(1, I).asArray(), 1);
        dgemv('No transpose', M - I, I, -ONE, A(I + 1, 1), LDA,
            X(1, I).asArray(), 1, ONE, X(I + 1, I).asArray(), 1);
        dgemv('No transpose', I - 1, N - I, ONE, A(1, I + 1), LDA,
            A(I, I + 1).asArray(), LDA, ZERO, X(1, I).asArray(), 1);
        dgemv('No transpose', M - I, I - 1, -ONE, X(I + 1, 1), LDX,
            X(1, I).asArray(), 1, ONE, X(I + 1, I).asArray(), 1);
        dscal(M - I, TAUP[I], X(I + 1, I).asArray(), 1);
      }
    }
  } else {
    // Reduce to lower bidiagonal form

    for (I = 1; I <= NB; I++) {
      // Update A(i,i:n)

      dgemv('No transpose', N - I + 1, I - 1, -ONE, Y(I, 1), LDY,
          A(I, 1).asArray(), LDA, ONE, A(I, I).asArray(), LDA);
      dgemv('Transpose', I - 1, N - I + 1, -ONE, A(1, I), LDA,
          X(I, 1).asArray(), LDX, ONE, A(I, I).asArray(), LDA);

      // Generate reflection P(i) to annihilate A(i,i+1:n)

      dlarfg(N - I + 1, A.box(I, I), A(I, min(I + 1, N)).asArray(), LDA,
          TAUP.box(I));
      D[I] = A[I][I];
      if (I < M) {
        A[I][I] = ONE;

        // Compute X(i+1:m,i)

        dgemv('No transpose', M - I, N - I + 1, ONE, A(I + 1, I), LDA,
            A(I, I).asArray(), LDA, ZERO, X(I + 1, I).asArray(), 1);
        dgemv('Transpose', N - I + 1, I - 1, ONE, Y(I, 1), LDY,
            A(I, I).asArray(), LDA, ZERO, X(1, I).asArray(), 1);
        dgemv('No transpose', M - I, I - 1, -ONE, A(I + 1, 1), LDA,
            X(1, I).asArray(), 1, ONE, X(I + 1, I).asArray(), 1);
        dgemv('No transpose', I - 1, N - I + 1, ONE, A(1, I), LDA,
            A(I, I).asArray(), LDA, ZERO, X(1, I).asArray(), 1);
        dgemv('No transpose', M - I, I - 1, -ONE, X(I + 1, 1), LDX,
            X(1, I).asArray(), 1, ONE, X(I + 1, I).asArray(), 1);
        dscal(M - I, TAUP[I], X(I + 1, I).asArray(), 1);

        // Update A(i+1:m,i)

        dgemv('No transpose', M - I, I - 1, -ONE, A(I + 1, 1), LDA,
            Y(I, 1).asArray(), LDY, ONE, A(I + 1, I).asArray(), 1);
        dgemv('No transpose', M - I, I, -ONE, X(I + 1, 1), LDX,
            A(1, I).asArray(), 1, ONE, A(I + 1, I).asArray(), 1);

        // Generate reflection Q(i) to annihilate A(i+2:m,i)

        dlarfg(M - I, A.box(I + 1, I), A(min(I + 2, M), I).asArray(), 1,
            TAUQ.box(I));
        E[I] = A[I + 1][I];
        A[I + 1][I] = ONE;

        // Compute Y(i+1:n,i)

        dgemv('Transpose', M - I, N - I, ONE, A(I + 1, I + 1), LDA,
            A(I + 1, I).asArray(), 1, ZERO, Y(I + 1, I).asArray(), 1);
        dgemv('Transpose', M - I, I - 1, ONE, A(I + 1, 1), LDA,
            A(I + 1, I).asArray(), 1, ZERO, Y(1, I).asArray(), 1);
        dgemv('No transpose', N - I, I - 1, -ONE, Y(I + 1, 1), LDY,
            Y(1, I).asArray(), 1, ONE, Y(I + 1, I).asArray(), 1);
        dgemv('Transpose', M - I, I, ONE, X(I + 1, 1), LDX,
            A(I + 1, I).asArray(), 1, ZERO, Y(1, I).asArray(), 1);
        dgemv('Transpose', I, N - I, -ONE, A(1, I + 1), LDA, Y(1, I).asArray(),
            1, ONE, Y(I + 1, I).asArray(), 1);
        dscal(N - I, TAUQ[I], Y(I + 1, I).asArray(), 1);
      }
    }
  }
}
