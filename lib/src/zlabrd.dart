import 'dart:math';

import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlacgv.dart';
import 'package:lapack/src/zlarfg.dart';

void zlabrd(
  final int M,
  final int N,
  final int NB,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> D_,
  final Array<double> E_,
  final Array<Complex> TAUQ_,
  final Array<Complex> TAUP_,
  final Matrix<Complex> X_,
  final int LDX,
  final Matrix<Complex> Y_,
  final int LDY,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final X = X_.having(ld: LDX);
  final Y = Y_.having(ld: LDY);
  final D = D_.having();
  final E = E_.having();
  final TAUQ = TAUQ_.having();
  final TAUP = TAUP_.having();
  int I;
  final ALPHA = Box(Complex.zero);

  // Quick return if possible

  if (M <= 0 || N <= 0) return;

  if (M >= N) {
    // Reduce to upper bidiagonal form

    for (I = 1; I <= NB; I++) {
      // Update A(i:m,i)

      zlacgv(I - 1, Y(I, 1).asArray(), LDY);
      zgemv('No transpose', M - I + 1, I - 1, -Complex.one, A(I, 1), LDA,
          Y(I, 1).asArray(), LDY, Complex.one, A(I, I).asArray(), 1);
      zlacgv(I - 1, Y(I, 1).asArray(), LDY);
      zgemv('No transpose', M - I + 1, I - 1, -Complex.one, X(I, 1), LDX,
          A(1, I).asArray(), 1, Complex.one, A(I, I).asArray(), 1);

      // Generate reflection Q(i) to annihilate A(i+1:m,i)

      ALPHA.value = A[I][I];
      zlarfg(M - I + 1, ALPHA, A(min(I + 1, M), I).asArray(), 1, TAUQ.box(I));
      D[I] = ALPHA.value.real;
      if (I < N) {
        A[I][I] = Complex.one;

        // Compute Y(i+1:n,i)

        zgemv('Conjugate transpose', M - I + 1, N - I, Complex.one, A(I, I + 1),
            LDA, A(I, I).asArray(), 1, Complex.zero, Y(I + 1, I).asArray(), 1);
        zgemv('Conjugate transpose', M - I + 1, I - 1, Complex.one, A(I, 1),
            LDA, A(I, I).asArray(), 1, Complex.zero, Y(1, I).asArray(), 1);
        zgemv('No transpose', N - I, I - 1, -Complex.one, Y(I + 1, 1), LDY,
            Y(1, I).asArray(), 1, Complex.one, Y(I + 1, I).asArray(), 1);
        zgemv('Conjugate transpose', M - I + 1, I - 1, Complex.one, X(I, 1),
            LDX, A(I, I).asArray(), 1, Complex.zero, Y(1, I).asArray(), 1);
        zgemv('Conjugate transpose', I - 1, N - I, -Complex.one, A(1, I + 1),
            LDA, Y(1, I).asArray(), 1, Complex.one, Y(I + 1, I).asArray(), 1);
        zscal(N - I, TAUQ[I], Y(I + 1, I).asArray(), 1);

        // Update A(i,i+1:n)

        zlacgv(N - I, A(I, I + 1).asArray(), LDA);
        zlacgv(I, A(I, 1).asArray(), LDA);
        zgemv('No transpose', N - I, I, -Complex.one, Y(I + 1, 1), LDY,
            A(I, 1).asArray(), LDA, Complex.one, A(I, I + 1).asArray(), LDA);
        zlacgv(I, A(I, 1).asArray(), LDA);
        zlacgv(I - 1, X(I, 1).asArray(), LDX);
        zgemv(
            'Conjugate transpose',
            I - 1,
            N - I,
            -Complex.one,
            A(1, I + 1),
            LDA,
            X(I, 1).asArray(),
            LDX,
            Complex.one,
            A(I, I + 1).asArray(),
            LDA);
        zlacgv(I - 1, X(I, 1).asArray(), LDX);

        // Generate reflection P(i) to annihilate A(i,i+2:n)

        ALPHA.value = A[I][I + 1];
        zlarfg(N - I, ALPHA, A(I, min(I + 2, N)).asArray(), LDA, TAUP.box(I));
        E[I] = ALPHA.value.real;
        A[I][I + 1] = Complex.one;

        // Compute X(i+1:m,i)

        zgemv('No transpose', M - I, N - I, Complex.one, A(I + 1, I + 1), LDA,
            A(I, I + 1).asArray(), LDA, Complex.zero, X(I + 1, I).asArray(), 1);
        zgemv('Conjugate transpose', N - I, I, Complex.one, Y(I + 1, 1), LDY,
            A(I, I + 1).asArray(), LDA, Complex.zero, X(1, I).asArray(), 1);
        zgemv('No transpose', M - I, I, -Complex.one, A(I + 1, 1), LDA,
            X(1, I).asArray(), 1, Complex.one, X(I + 1, I).asArray(), 1);
        zgemv('No transpose', I - 1, N - I, Complex.one, A(1, I + 1), LDA,
            A(I, I + 1).asArray(), LDA, Complex.zero, X(1, I).asArray(), 1);
        zgemv('No transpose', M - I, I - 1, -Complex.one, X(I + 1, 1), LDX,
            X(1, I).asArray(), 1, Complex.one, X(I + 1, I).asArray(), 1);
        zscal(M - I, TAUP[I], X(I + 1, I).asArray(), 1);
        zlacgv(N - I, A(I, I + 1).asArray(), LDA);
      }
    }
  } else {
    // Reduce to lower bidiagonal form

    for (I = 1; I <= NB; I++) {
      // Update A(i,i:n)

      zlacgv(N - I + 1, A(I, I).asArray(), LDA);
      zlacgv(I - 1, A(I, 1).asArray(), LDA);
      zgemv('No transpose', N - I + 1, I - 1, -Complex.one, Y(I, 1), LDY,
          A(I, 1).asArray(), LDA, Complex.one, A(I, I).asArray(), LDA);
      zlacgv(I - 1, A(I, 1).asArray(), LDA);
      zlacgv(I - 1, X(I, 1).asArray(), LDX);
      zgemv('Conjugate transpose', I - 1, N - I + 1, -Complex.one, A(1, I), LDA,
          X(I, 1).asArray(), LDX, Complex.one, A(I, I).asArray(), LDA);
      zlacgv(I - 1, X(I, 1).asArray(), LDX);

      // Generate reflection P(i) to annihilate A(i,i+1:n)

      ALPHA.value = A[I][I];
      zlarfg(N - I + 1, ALPHA, A(I, min(I + 1, N)).asArray(), LDA, TAUP.box(I));
      D[I] = ALPHA.value.real;
      if (I < M) {
        A[I][I] = Complex.one;

        // Compute X(i+1:m,i)

        zgemv('No transpose', M - I, N - I + 1, Complex.one, A(I + 1, I), LDA,
            A(I, I).asArray(), LDA, Complex.zero, X(I + 1, I).asArray(), 1);
        zgemv('Conjugate transpose', N - I + 1, I - 1, Complex.one, Y(I, 1),
            LDY, A(I, I).asArray(), LDA, Complex.zero, X(1, I).asArray(), 1);
        zgemv('No transpose', M - I, I - 1, -Complex.one, A(I + 1, 1), LDA,
            X(1, I).asArray(), 1, Complex.one, X(I + 1, I).asArray(), 1);
        zgemv('No transpose', I - 1, N - I + 1, Complex.one, A(1, I), LDA,
            A(I, I).asArray(), LDA, Complex.zero, X(1, I).asArray(), 1);
        zgemv('No transpose', M - I, I - 1, -Complex.one, X(I + 1, 1), LDX,
            X(1, I).asArray(), 1, Complex.one, X(I + 1, I).asArray(), 1);
        zscal(M - I, TAUP[I], X(I + 1, I).asArray(), 1);
        zlacgv(N - I + 1, A(I, I).asArray(), LDA);

        // Update A(i+1:m,i)

        zlacgv(I - 1, Y(I, 1).asArray(), LDY);
        zgemv('No transpose', M - I, I - 1, -Complex.one, A(I + 1, 1), LDA,
            Y(I, 1).asArray(), LDY, Complex.one, A(I + 1, I).asArray(), 1);
        zlacgv(I - 1, Y(I, 1).asArray(), LDY);
        zgemv('No transpose', M - I, I, -Complex.one, X(I + 1, 1), LDX,
            A(1, I).asArray(), 1, Complex.one, A(I + 1, I).asArray(), 1);

        // Generate reflection Q(i) to annihilate A(i+2:m,i)

        ALPHA.value = A[I + 1][I];
        zlarfg(M - I, ALPHA, A(min(I + 2, M), I).asArray(), 1, TAUQ.box(I));
        E[I] = ALPHA.value.real;
        A[I + 1][I] = Complex.one;

        // Compute Y(i+1:n,i)

        zgemv(
            'Conjugate transpose',
            M - I,
            N - I,
            Complex.one,
            A(I + 1, I + 1),
            LDA,
            A(I + 1, I).asArray(),
            1,
            Complex.zero,
            Y(I + 1, I).asArray(),
            1);
        zgemv('Conjugate transpose', M - I, I - 1, Complex.one, A(I + 1, 1),
            LDA, A(I + 1, I).asArray(), 1, Complex.zero, Y(1, I).asArray(), 1);
        zgemv('No transpose', N - I, I - 1, -Complex.one, Y(I + 1, 1), LDY,
            Y(1, I).asArray(), 1, Complex.one, Y(I + 1, I).asArray(), 1);
        zgemv('Conjugate transpose', M - I, I, Complex.one, X(I + 1, 1), LDX,
            A(I + 1, I).asArray(), 1, Complex.zero, Y(1, I).asArray(), 1);
        zgemv('Conjugate transpose', I, N - I, -Complex.one, A(1, I + 1), LDA,
            Y(1, I).asArray(), 1, Complex.one, Y(I + 1, I).asArray(), 1);
        zscal(N - I, TAUQ[I], Y(I + 1, I).asArray(), 1);
      } else {
        zlacgv(N - I + 1, A(I, I).asArray(), LDA);
      }
    }
  }
}
