import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zdotc.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/blas/zhemv.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlacgv.dart';
import 'package:lapack/src/zlarfg.dart';

void zlatrd(
  final String UPLO,
  final int N,
  final int NB,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> E_,
  final Array<Complex> TAU_,
  final Matrix<Complex> W_,
  final int LDW,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final E = E_.having();
  final TAU = TAU_.having();
  final W = W_.having(ld: LDW);
  const HALF = Complex(0.5, 0.0);
  int I, IW;
  final ALPHA = Box(Complex.zero);

  // Quick return if possible

  if (N <= 0) return;

  if (lsame(UPLO, 'U')) {
    // Reduce last NB columns of upper triangle

    for (I = N; I >= N - NB + 1; I--) {
      IW = I - N + NB;
      if (I < N) {
        // Update A(1:i,i)

        A[I][I] = A[I][I].real.toComplex();
        zlacgv(N - I, W(I, IW + 1).asArray(), LDW);
        zgemv('No transpose', I, N - I, -Complex.one, A(1, I + 1), LDA,
            W(I, IW + 1).asArray(), LDW, Complex.one, A(1, I).asArray(), 1);
        zlacgv(N - I, W(I, IW + 1).asArray(), LDW);
        zlacgv(N - I, A(I, I + 1).asArray(), LDA);
        zgemv('No transpose', I, N - I, -Complex.one, W(1, IW + 1), LDW,
            A(I, I + 1).asArray(), LDA, Complex.one, A(1, I).asArray(), 1);
        zlacgv(N - I, A(I, I + 1).asArray(), LDA);
        A[I][I] = A[I][I].real.toComplex();
      }
      if (I > 1) {
        // Generate elementary reflector H(i) to annihilate
        // A(1:i-2,i)

        ALPHA.value = A[I - 1][I];
        zlarfg(I - 1, ALPHA, A(1, I).asArray(), 1, TAU(I - 1));
        E[I - 1] = ALPHA.value.real;
        A[I - 1][I] = Complex.one;

        // Compute W(1:i-1,i)

        zhemv('Upper', I - 1, Complex.one, A, LDA, A(1, I).asArray(), 1,
            Complex.zero, W(1, IW).asArray(), 1);
        if (I < N) {
          zgemv(
              'Conjugate transpose',
              I - 1,
              N - I,
              Complex.one,
              W(1, IW + 1),
              LDW,
              A(1, I).asArray(),
              1,
              Complex.zero,
              W(I + 1, IW).asArray(),
              1);
          zgemv('No transpose', I - 1, N - I, -Complex.one, A(1, I + 1), LDA,
              W(I + 1, IW).asArray(), 1, Complex.one, W(1, IW).asArray(), 1);
          zgemv(
              'Conjugate transpose',
              I - 1,
              N - I,
              Complex.one,
              A(1, I + 1),
              LDA,
              A(1, I).asArray(),
              1,
              Complex.zero,
              W(I + 1, IW).asArray(),
              1);
          zgemv('No transpose', I - 1, N - I, -Complex.one, W(1, IW + 1), LDW,
              W(I + 1, IW).asArray(), 1, Complex.one, W(1, IW).asArray(), 1);
        }
        zscal(I - 1, TAU[I - 1], W(1, IW).asArray(), 1);
        ALPHA.value = -HALF *
            TAU[I - 1] *
            zdotc(I - 1, W(1, IW).asArray(), 1, A(1, I).asArray(), 1);
        zaxpy(I - 1, ALPHA.value, A(1, I).asArray(), 1, W(1, IW).asArray(), 1);
      }
    }
  } else {
    // Reduce first NB columns of lower triangle

    for (I = 1; I <= NB; I++) {
      // Update A(i:n,i)

      A[I][I] = A[I][I].real.toComplex();
      zlacgv(I - 1, W(I, 1).asArray(), LDW);
      zgemv('No transpose', N - I + 1, I - 1, -Complex.one, A(I, 1), LDA,
          W(I, 1).asArray(), LDW, Complex.one, A(I, I).asArray(), 1);
      zlacgv(I - 1, W(I, 1).asArray(), LDW);
      zlacgv(I - 1, A(I, 1).asArray(), LDA);
      zgemv('No transpose', N - I + 1, I - 1, -Complex.one, W(I, 1), LDW,
          A(I, 1).asArray(), LDA, Complex.one, A(I, I).asArray(), 1);
      zlacgv(I - 1, A(I, 1).asArray(), LDA);
      A[I][I] = A[I][I].real.toComplex();
      if (I < N) {
        // Generate elementary reflector H(i) to annihilate
        // A(i+2:n,i)

        ALPHA.value = A[I + 1][I];
        zlarfg(N - I, ALPHA, A(min(I + 2, N), I).asArray(), 1, TAU(I));
        E[I] = ALPHA.value.real;
        A[I + 1][I] = Complex.one;

        // Compute W(i+1:n,i)

        zhemv('Lower', N - I, Complex.one, A(I + 1, I + 1), LDA,
            A(I + 1, I).asArray(), 1, Complex.zero, W(I + 1, I).asArray(), 1);
        zgemv('Conjugate transpose', N - I, I - 1, Complex.one, W(I + 1, 1),
            LDW, A(I + 1, I).asArray(), 1, Complex.zero, W(1, I).asArray(), 1);
        zgemv('No transpose', N - I, I - 1, -Complex.one, A(I + 1, 1), LDA,
            W(1, I).asArray(), 1, Complex.one, W(I + 1, I).asArray(), 1);
        zgemv('Conjugate transpose', N - I, I - 1, Complex.one, A(I + 1, 1),
            LDA, A(I + 1, I).asArray(), 1, Complex.zero, W(1, I).asArray(), 1);
        zgemv('No transpose', N - I, I - 1, -Complex.one, W(I + 1, 1), LDW,
            W(1, I).asArray(), 1, Complex.one, W(I + 1, I).asArray(), 1);
        zscal(N - I, TAU[I], W(I + 1, I).asArray(), 1);
        ALPHA.value = -HALF *
            TAU[I] *
            zdotc(N - I, W(I + 1, I).asArray(), 1, A(I + 1, I).asArray(), 1);
        zaxpy(N - I, ALPHA.value, A(I + 1, I).asArray(), 1,
            W(I + 1, I).asArray(), 1);
      }
    }
  }
}
