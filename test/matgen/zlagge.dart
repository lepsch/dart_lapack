// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dznrm2.dart';
import 'package:dart_lapack/src/blas/zgemv.dart';
import 'package:dart_lapack/src/blas/zgerc.dart';
import 'package:dart_lapack/src/blas/zscal.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlacgv.dart';
import 'package:dart_lapack/src/zlarnv.dart';

void zlagge(
  final int M,
  final int N,
  final int KL,
  final int KU,
  final Array<double> D_,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> ISEED_,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
  final D = D_.having();
  final A = A_.having(ld: LDA);
  final ISEED = ISEED_.having(length: 4);
  final WORK = WORK_.having();
  const ZERO = 0.0;
  int I, J;
  double WN;
  Complex TAU, WA, WB;

  // Test the input arguments

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (KL < 0 || KL > M - 1) {
    INFO.value = -3;
  } else if (KU < 0 || KU > N - 1) {
    INFO.value = -4;
  } else if (LDA < max(1, M)) {
    INFO.value = -7;
  }
  if (INFO.value < 0) {
    xerbla('ZLAGGE', -INFO.value);
    return;
  }

  // initialize A to diagonal matrix

  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      A[I][J] = Complex.zero;
    }
  }
  for (I = 1; I <= min(M, N); I++) {
    A[I][I] = D[I].toComplex();
  }

  // Quick exit if the user wants a diagonal matrix

  if ((KL == 0) && (KU == 0)) return;

  // pre- and post-multiply A by random unitary matrices

  for (I = min(M, N); I >= 1; I--) {
    if (I < M) {
      // generate random reflection

      zlarnv(3, ISEED, M - I + 1, WORK);
      WN = dznrm2(M - I + 1, WORK, 1);
      WA = (WN / WORK[1].abs()).toComplex() * WORK[1];
      if (WN == ZERO) {
        TAU = Complex.zero;
      } else {
        WB = WORK[1] + WA;
        zscal(M - I, Complex.one / WB, WORK(2), 1);
        WORK[1] = Complex.one;
        TAU = (WB / WA).real.toComplex();
      }

      // multiply A(i:m,i:n) by random reflection from the left

      zgemv('Conjugate transpose', M - I + 1, N - I + 1, Complex.one, A(I, I),
          LDA, WORK, 1, Complex.zero, WORK(M + 1), 1);
      zgerc(M - I + 1, N - I + 1, -TAU, WORK, 1, WORK(M + 1), 1, A(I, I), LDA);
    }
    if (I < N) {
      // generate random reflection

      zlarnv(3, ISEED, N - I + 1, WORK);
      WN = dznrm2(N - I + 1, WORK, 1);
      WA = (WN / WORK[1].abs()).toComplex() * WORK[1];
      if (WN == ZERO) {
        TAU = Complex.zero;
      } else {
        WB = WORK[1] + WA;
        zscal(N - I, Complex.one / WB, WORK(2), 1);
        WORK[1] = Complex.one;
        TAU = (WB / WA).real.toComplex();
      }

      // multiply A(i:m,i:n) by random reflection from the right

      zgemv('No transpose', M - I + 1, N - I + 1, Complex.one, A(I, I), LDA,
          WORK, 1, Complex.zero, WORK(N + 1), 1);
      zgerc(M - I + 1, N - I + 1, -TAU, WORK(N + 1), 1, WORK, 1, A(I, I), LDA);
    }
  }

  // Reduce number of subdiagonals to KL and number of superdiagonals
  // to KU

  for (I = 1; I <= max(M - 1 - KL, N - 1 - KU); I++) {
    if (KL <= KU) {
      // annihilate subdiagonal elements first (necessary if KL = 0)

      if (I <= min(M - 1 - KL, N)) {
        // generate reflection to annihilate A(kl+i+1:m,i)

        WN = dznrm2(M - KL - I + 1, A(KL + I, I).asArray(), 1);
        WA = (WN / A[KL + I][I].abs()).toComplex() * A[KL + I][I];
        if (WN == ZERO) {
          TAU = Complex.zero;
        } else {
          WB = A[KL + I][I] + WA;
          zscal(M - KL - I, Complex.one / WB, A(KL + I + 1, I).asArray(), 1);
          A[KL + I][I] = Complex.one;
          TAU = (WB / WA).real.toComplex();
        }

        // apply reflection to A(kl+i:m,i+1:n) from the left

        zgemv(
            'Conjugate transpose',
            M - KL - I + 1,
            N - I,
            Complex.one,
            A(KL + I, I + 1),
            LDA,
            A(KL + I, I).asArray(),
            1,
            Complex.zero,
            WORK,
            1);
        zgerc(M - KL - I + 1, N - I, -TAU, A(KL + I, I).asArray(), 1, WORK, 1,
            A(KL + I, I + 1), LDA);
        A[KL + I][I] = -WA;
      }

      if (I <= min(N - 1 - KU, M)) {
        // generate reflection to annihilate A(i,ku+i+1:n)

        WN = dznrm2(N - KU - I + 1, A(I, KU + I).asArray(), LDA);
        WA = (WN / A[I][KU + I].abs()).toComplex() * A[I][KU + I];
        if (WN == ZERO) {
          TAU = Complex.zero;
        } else {
          WB = A[I][KU + I] + WA;
          zscal(N - KU - I, Complex.one / WB, A(I, KU + I + 1).asArray(), LDA);
          A[I][KU + I] = Complex.one;
          TAU = (WB / WA).real.toComplex();
        }

        // apply reflection to A(i+1:m,ku+i:n) from the right

        zlacgv(N - KU - I + 1, A(I, KU + I).asArray(), LDA);
        zgemv(
            'No transpose',
            M - I,
            N - KU - I + 1,
            Complex.one,
            A(I + 1, KU + I),
            LDA,
            A(I, KU + I).asArray(),
            LDA,
            Complex.zero,
            WORK,
            1);
        zgerc(M - I, N - KU - I + 1, -TAU, WORK, 1, A(I, KU + I).asArray(), LDA,
            A(I + 1, KU + I), LDA);
        A[I][KU + I] = -WA;
      }
    } else {
      // annihilate superdiagonal elements first (necessary if
      // KU = 0)

      if (I <= min(N - 1 - KU, M)) {
        // generate reflection to annihilate A(i,ku+i+1:n)

        WN = dznrm2(N - KU - I + 1, A(I, KU + I).asArray(), LDA);
        WA = (WN / A[I][KU + I].abs()).toComplex() * A[I][KU + I];
        if (WN == ZERO) {
          TAU = Complex.zero;
        } else {
          WB = A[I][KU + I] + WA;
          zscal(N - KU - I, Complex.one / WB, A(I, KU + I + 1).asArray(), LDA);
          A[I][KU + I] = Complex.one;
          TAU = (WB / WA).real.toComplex();
        }

        // apply reflection to A(i+1:m,ku+i:n) from the right

        zlacgv(N - KU - I + 1, A(I, KU + I).asArray(), LDA);
        zgemv(
            'No transpose',
            M - I,
            N - KU - I + 1,
            Complex.one,
            A(I + 1, KU + I),
            LDA,
            A(I, KU + I).asArray(),
            LDA,
            Complex.zero,
            WORK,
            1);
        zgerc(M - I, N - KU - I + 1, -TAU, WORK, 1, A(I, KU + I).asArray(), LDA,
            A(I + 1, KU + I), LDA);
        A[I][KU + I] = -WA;
      }

      if (I <= min(M - 1 - KL, N)) {
        // generate reflection to annihilate A(kl+i+1:m,i)

        WN = dznrm2(M - KL - I + 1, A(KL + I, I).asArray(), 1);
        WA = (WN / A[KL + I][I].abs()).toComplex() * A[KL + I][I];
        if (WN == ZERO) {
          TAU = Complex.zero;
        } else {
          WB = A[KL + I][I] + WA;
          zscal(M - KL - I, Complex.one / WB, A(KL + I + 1, I).asArray(), 1);
          A[KL + I][I] = Complex.one;
          TAU = (WB / WA).real.toComplex();
        }

        // apply reflection to A(kl+i:m,i+1:n) from the left

        zgemv(
            'Conjugate transpose',
            M - KL - I + 1,
            N - I,
            Complex.one,
            A(KL + I, I + 1),
            LDA,
            A(KL + I, I).asArray(),
            1,
            Complex.zero,
            WORK,
            1);
        zgerc(M - KL - I + 1, N - I, -TAU, A(KL + I, I).asArray(), 1, WORK, 1,
            A(KL + I, I + 1), LDA);
        A[KL + I][I] = -WA;
      }
    }

    if (I <= N) {
      for (J = KL + I + 1; J <= M; J++) {
        A[J][I] = Complex.zero;
      }
    }

    if (I <= M) {
      for (J = KU + I + 1; J <= N; J++) {
        A[I][J] = Complex.zero;
      }
    }
  }
}
