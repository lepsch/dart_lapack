// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zaxpy.dart';
import 'package:dart_lapack/src/blas/zdotc.dart';
import 'package:dart_lapack/src/blas/zhpmv.dart';
import 'package:dart_lapack/src/blas/zhpr2.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlarfg.dart';

void zhptrd(
  final String UPLO,
  final int N,
  final Array<Complex> AP_,
  final Array<double> D_,
  final Array<double> E_,
  final Array<Complex> TAU_,
  final Box<int> INFO,
) {
  final AP = AP_.having();
  final TAU = TAU_.having();
  final D = D_.having();
  final E = E_.having();
  const HALF = Complex(0.5, 0.0);
  bool UPPER;
  int I = 0, I1, I1I1, II;
  final ALPHA = Box(Complex.zero), TAUI = Box(Complex.zero);

  // Test the input parameters

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  }
  if (INFO.value != 0) {
    xerbla('ZHPTRD', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N <= 0) return;

  if (UPPER) {
    // Reduce the upper triangle of A.
    // I1 is the index in AP of A(1,I+1).

    I1 = N * (N - 1) ~/ 2 + 1;
    AP[I1 + N - 1] = AP[I1 + N - 1].real.toComplex();
    for (I = N - 1; I >= 1; I--) {
      // Generate elementary reflector H(i) = I - tau * v * v**H
      // to annihilate A(1:i-1,i+1)

      ALPHA.value = AP[I1 + I - 1];
      zlarfg(I, ALPHA, AP(I1), 1, TAUI);
      E[I] = ALPHA.value.real;

      if (TAUI.value != Complex.zero) {
        // Apply H(i) from both sides to A(1:i,1:i)

        AP[I1 + I - 1] = Complex.one;

        // Compute  y := tau * A * v  storing y in TAU(1:i)

        zhpmv(UPLO, I, TAUI.value, AP, AP(I1), 1, Complex.zero, TAU, 1);

        // Compute  w := y - 1/2 * tau * (y**H *v) * v

        ALPHA.value = -HALF * TAUI.value * zdotc(I, TAU, 1, AP(I1), 1);
        zaxpy(I, ALPHA.value, AP(I1), 1, TAU, 1);

        // Apply the transformation as a rank-2 update:
        //    A := A - v * w**H - w * v**H

        zhpr2(UPLO, I, -Complex.one, AP(I1), 1, TAU, 1, AP);
      }
      AP[I1 + I - 1] = E[I].toComplex();
      D[I + 1] = AP[I1 + I].real;
      TAU[I] = TAUI.value;
      I1 -= I;
    }
    D[1] = AP[1].real;
  } else {
    // Reduce the lower triangle of A. II is the index in AP of
    // A(i,i) and I1I1 is the index of A(i+1,i+1).

    II = 1;
    AP[1] = AP[1].real.toComplex();
    for (I = 1; I <= N - 1; I++) {
      I1I1 = II + N - I + 1;

      // Generate elementary reflector H(i) = I - tau * v * v**H
      // to annihilate A(i+2:n,i)

      ALPHA.value = AP[II + 1];
      zlarfg(N - I, ALPHA, AP(II + 2), 1, TAUI);
      E[I] = ALPHA.value.real;

      if (TAUI.value != Complex.zero) {
        // Apply H(i) from both sides to A(i+1:n,i+1:n)

        AP[II + 1] = Complex.one;

        // Compute  y := tau * A * v  storing y in TAU(i:n-1)

        zhpmv(UPLO, N - I, TAUI.value, AP(I1I1), AP(II + 1), 1, Complex.zero,
            TAU(I), 1);

        // Compute  w := y - 1/2 * tau * (y**H *v) * v

        ALPHA.value =
            -HALF * TAUI.value * zdotc(N - I, TAU(I), 1, AP(II + 1), 1);
        zaxpy(N - I, ALPHA.value, AP(II + 1), 1, TAU(I), 1);

        // Apply the transformation as a rank-2 update:
        //    A := A - v * w**H - w * v**H

        zhpr2(UPLO, N - I, -Complex.one, AP(II + 1), 1, TAU(I), 1, AP(I1I1));
      }
      AP[II + 1] = E[I].toComplex();
      D[I] = AP[II].real;
      TAU[I] = TAUI.value;
      II = I1I1;
    }
    D[N] = AP[II].real;
  }
}
