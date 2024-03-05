import 'dart:math';

import 'package:lapack/src/blas/dznrm2.dart';
import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zdotc.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/blas/zgerc.dart';
import 'package:lapack/src/blas/zhemv.dart';
import 'package:lapack/src/blas/zher2.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlarnv.dart';

void zlaghe(
  final int N,
  final int K,
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

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const HALF = Complex(0.5, 0.0);
  const ZERO = 0.0;
  int I, J;
  double WN;
  Complex ALPHA, TAU, WA, WB;

  // Test the input arguments

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
  } else if (K < 0 || K > N - 1) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  }
  if (INFO.value < 0) {
    xerbla('ZLAGHE', -INFO.value);
    return;
  }

  // initialize lower triangle of A to diagonal matrix

  for (J = 1; J <= N; J++) {
    // 20
    for (I = J + 1; I <= N; I++) {
      // 10
      A[I][J] = Complex.zero;
    } // 10
  } // 20
  for (I = 1; I <= N; I++) {
    // 30
    A[I][I] = D[I].toComplex();
  } // 30

  // Generate lower triangle of hermitian matrix

  for (I = N - 1; I >= 1; I--) {
    // 40

    // generate random reflection

    zlarnv(3, ISEED, N - I + 1, WORK);
    WN = dznrm2(N - I + 1, WORK, 1);
    WA = (WN / (WORK[1]).abs()).toComplex() * WORK[1];
    if (WN == ZERO) {
      TAU = Complex.zero;
    } else {
      WB = WORK[1] + WA;
      zscal(N - I, Complex.one / WB, WORK(2), 1);
      WORK[1] = Complex.one;
      TAU = (WB / WA).toDouble().toComplex();
    }

    // apply random reflection to A(i:n,i:n) from the left
    // and the right

    // compute  y := tau * A * u

    zhemv('Lower', N - I + 1, TAU, A(I, I), LDA, WORK, 1, Complex.zero,
        WORK(N + 1), 1);

    // compute  v := y - 1/2 * tau * ( y, u ) * u

    ALPHA = -HALF * TAU * zdotc(N - I + 1, WORK(N + 1), 1, WORK, 1);
    zaxpy(N - I + 1, ALPHA, WORK, 1, WORK(N + 1), 1);

    // apply the transformation as a rank-2 update to A(i:n,i:n)

    zher2('Lower', N - I + 1, -Complex.one, WORK, 1, WORK(N + 1), 1, A(I, I),
        LDA);
  } // 40

  // Reduce number of subdiagonals to K

  for (I = 1; I <= N - 1 - K; I++) {
    // 60

    // generate reflection to annihilate A(k+i+1:n,i)

    WN = dznrm2(N - K - I + 1, A(K + I, I).asArray(), 1);
    WA = (WN / (A[K + I][I]).abs()).toComplex() * A[K + I][I];
    if (WN == ZERO) {
      TAU = Complex.zero;
    } else {
      WB = A[K + I][I] + WA;
      zscal(N - K - I, Complex.one / WB, A(K + I + 1, I).asArray(), 1);
      A[K + I][I] = Complex.one;
      TAU = (WB / WA).toDouble().toComplex();
    }

    // apply reflection to A(k+i:n,i+1:k+i-1) from the left

    zgemv('Conjugate transpose', N - K - I + 1, K - 1, Complex.one,
        A(K + I, I + 1), LDA, A(K + I, I).asArray(), 1, Complex.zero, WORK, 1);
    zgerc(N - K - I + 1, K - 1, -TAU, A(K + I, I).asArray(), 1, WORK, 1,
        A(K + I, I + 1), LDA);

    // apply reflection to A(k+i:n,k+i:n) from the left and the right

    // compute  y := tau * A * u

    zhemv('Lower', N - K - I + 1, TAU, A(K + I, K + I), LDA,
        A(K + I, I).asArray(), 1, Complex.zero, WORK, 1);

    // compute  v := y - 1/2 * tau * ( y, u ) * u

    ALPHA =
        -HALF * TAU * zdotc(N - K - I + 1, WORK, 1, A(K + I, I).asArray(), 1);
    zaxpy(N - K - I + 1, ALPHA, A(K + I, I).asArray(), 1, WORK, 1);

    // apply hermitian rank-2 update to A(k+i:n,k+i:n)

    zher2('Lower', N - K - I + 1, -Complex.one, A(K + I, I).asArray(), 1, WORK,
        1, A(K + I, K + I), LDA);

    A[K + I][I] = -WA;
    for (J = K + I + 1; J <= N; J++) {
      // 50
      A[J][I] = Complex.zero;
    } // 50
  } // 60

  // Store full hermitian matrix

  for (J = 1; J <= N; J++) {
    // 80
    for (I = J + 1; I <= N; I++) {
      // 70
      A[J][I] = A[I][J].abs().toComplex();
    } // 70
  } // 80
}
