import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacgv.dart';
import 'package:lapack/src/zlarf.dart';
import 'package:lapack/src/zlarfg.dart';

void zgebd2(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> D_,
  final Array<double> E_,
  final Array<Complex> TAUQ_,
  final Array<Complex> TAUP_,
  final Array<Complex> WORK_,
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
  int I;
  final ALPHA = Box(Complex.zero);

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
    xerbla('ZGEBD2', -INFO.value);
    return;
  }

  if (M >= N) {
    // Reduce to upper bidiagonal form

    for (I = 1; I <= N; I++) {
      // Generate elementary reflector H(i) to annihilate A(i+1:m,i)

      ALPHA.value = A[I][I];
      zlarfg(M - I + 1, ALPHA, A(min(I + 1, M), I).asArray(), 1, TAUQ.box(I));
      D[I] = ALPHA.value.real;
      A[I][I] = Complex.one;

      // Apply H(i)**H to A(i:m,i+1:n) from the left

      if (I < N) {
        zlarf('Left', M - I + 1, N - I, A(I, I).asArray(), 1,
            TAUQ[I].conjugate(), A(I, I + 1), LDA, WORK);
      }
      A[I][I] = D[I].toComplex();

      if (I < N) {
        // Generate elementary reflector G(i) to annihilate
        // A(i,i+2:n)

        zlacgv(N - I, A(I, I + 1).asArray(), LDA);
        ALPHA.value = A[I][I + 1];
        zlarfg(N - I, ALPHA, A(I, min(I + 2, N)).asArray(), LDA, TAUP.box(I));
        E[I] = ALPHA.value.real;
        A[I][I + 1] = Complex.one;

        // Apply G(i) to A(i+1:m,i+1:n) from the right

        zlarf('Right', M - I, N - I, A(I, I + 1).asArray(), LDA, TAUP[I],
            A(I + 1, I + 1), LDA, WORK);
        zlacgv(N - I, A(I, I + 1).asArray(), LDA);
        A[I][I + 1] = E[I].toComplex();
      } else {
        TAUP[I] = Complex.zero;
      }
    }
  } else {
    // Reduce to lower bidiagonal form

    for (I = 1; I <= M; I++) {
      // Generate elementary reflector G(i) to annihilate A(i,i+1:n)

      zlacgv(N - I + 1, A(I, I).asArray(), LDA);
      ALPHA.value = A[I][I];
      zlarfg(N - I + 1, ALPHA, A(I, min(I + 1, N)).asArray(), LDA, TAUP.box(I));
      D[I] = ALPHA.value.real;
      A[I][I] = Complex.one;

      // Apply G(i) to A(i+1:m,i:n) from the right

      if (I < M) {
        zlarf('Right', M - I, N - I + 1, A(I, I).asArray(), LDA, TAUP[I],
            A(I + 1, I), LDA, WORK);
      }
      zlacgv(N - I + 1, A(I, I).asArray(), LDA);
      A[I][I] = D[I].toComplex();

      if (I < M) {
        // Generate elementary reflector H(i) to annihilate
        // A(i+2:m,i)

        ALPHA.value = A[I + 1][I];
        zlarfg(M - I, ALPHA, A(min(I + 2, M), I).asArray(), 1, TAUQ.box(I));
        E[I] = ALPHA.value.real;
        A[I + 1][I] = Complex.one;

        // Apply H(i)**H to A(i+1:m,i+1:n) from the left

        zlarf('Left', M - I, N - I, A(I + 1, I).asArray(), 1,
            TAUQ[I].conjugate(), A(I + 1, I + 1), LDA, WORK);
        A[I + 1][I] = E[I].toComplex();
      } else {
        TAUQ[I] = Complex.zero;
      }
    }
  }
}
