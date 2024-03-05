import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zdotc.dart';
import 'package:lapack/src/blas/zhemv.dart';
import 'package:lapack/src/blas/zher2.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlarfg.dart';

void zhetd2(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> D_,
  final Array<double> E_,
  final Array<Complex> TAU_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final D = D_.having();
  final E = E_.having();
  const HALF = Complex(0.5, 0.0);
  bool UPPER;
  int I;
  final ALPHA = Box(Complex.zero), TAUI = Box(Complex.zero);

  // Test the input parameters

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('ZHETD2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N <= 0) return;

  if (UPPER) {
    // Reduce the upper triangle of A

    A[N][N] = A[N][N].real.toComplex();
    for (I = N - 1; I >= 1; I--) {
      // 10

      // Generate elementary reflector H(i) = I - tau * v * v**H
      // to annihilate A(1:i-1,i+1)

      ALPHA.value = A[I][I + 1];
      zlarfg(I, ALPHA, A(1, I + 1).asArray(), 1, TAUI);
      E[I] = ALPHA.value.toDouble();

      if (TAUI.value != Complex.zero) {
        // Apply H(i) from both sides to A(1:i,1:i)

        A[I][I + 1] = Complex.one;

        // Compute  x := tau * A * v  storing x in TAU(1:i)

        zhemv(UPLO, I, TAUI.value, A, LDA, A(1, I + 1).asArray(), 1,
            Complex.zero, TAU, 1);

        // Compute  w := x - 1/2 * tau * (x**H * v) * v

        ALPHA.value =
            -HALF * TAUI.value * zdotc(I, TAU, 1, A(1, I + 1).asArray(), 1);
        zaxpy(I, ALPHA.value, A(1, I + 1).asArray(), 1, TAU, 1);

        // Apply the transformation as a rank-2 update:
        //    A := A - v * w**H - w * v**H

        zher2(UPLO, I, -Complex.one, A(1, I + 1).asArray(), 1, TAU, 1, A, LDA);
      } else {
        A[I][I] = A[I][I].real.toComplex();
      }
      A[I][I + 1] = E[I].toComplex();
      D[I + 1] = A[I + 1][I + 1].toDouble();
      TAU[I] = TAUI.value;
    } // 10
    D[1] = A[1][1].toDouble();
  } else {
    // Reduce the lower triangle of A

    A[1][1] = A[1][1].real.toComplex();
    for (I = 1; I <= N - 1; I++) {
      // 20

      // Generate elementary reflector H(i) = I - tau * v * v**H
      // to annihilate A(i+2:n,i)

      ALPHA.value = A[I + 1][I];
      zlarfg(N - I, ALPHA, A(min(I + 2, N), I).asArray(), 1, TAUI);
      E[I] = ALPHA.value.toDouble();

      if (TAUI.value != Complex.zero) {
        // Apply H(i) from both sides to A(i+1:n,i+1:n)

        A[I + 1][I] = Complex.one;

        // Compute  x := tau * A * v  storing y in TAU(i:n-1)

        zhemv(UPLO, N - I, TAUI.value, A(I + 1, I + 1), LDA,
            A(I + 1, I).asArray(), 1, Complex.zero, TAU(I), 1);

        // Compute  w := x - 1/2 * tau * (x**H * v) * v

        ALPHA.value = -HALF *
            TAUI.value *
            zdotc(N - I, TAU(I), 1, A(I + 1, I).asArray(), 1);
        zaxpy(N - I, ALPHA.value, A(I + 1, I).asArray(), 1, TAU(I), 1);

        // Apply the transformation as a rank-2 update:
        //    A := A - v * w**H - w * v**H

        zher2(UPLO, N - I, -Complex.one, A(I + 1, I).asArray(), 1, TAU(I), 1,
            A(I + 1, I + 1), LDA);
      } else {
        A[I + 1][I + 1] = A[I + 1][I + 1].real.toComplex();
      }
      A[I + 1][I] = E[I].toComplex();
      D[I] = A[I][I].toDouble();
      TAU[I] = TAUI.value;
    } // 20
    D[N] = (A[N][N]).toDouble();
  }
}
