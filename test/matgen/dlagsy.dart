import 'dart:math';

import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dger.dart';
import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dsymv.dart';
import 'package:lapack/src/blas/dsyr2.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlarnv.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlagsy(
  final int N,
  final int K,
  final Array<double> D_,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> ISEED_,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.dim();
  final A = A_.dim(LDA);
  final ISEED = ISEED_.dim();
  final WORK = WORK_.dim();
  const ZERO = 0.0, ONE = 1.0, HALF = 0.5;
  int I, J;
  double ALPHA, TAU, WA, WB, WN;

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
    xerbla('DLAGSY', -INFO.value);
    return;
  }

  // initialize lower triangle of A to diagonal matrix

  for (J = 1; J <= N; J++) {
    for (I = J + 1; I <= N; I++) {
      A[I][J] = ZERO;
    }
  }
  for (I = 1; I <= N; I++) {
    A[I][I] = D[I];
  }

  // Generate lower triangle of symmetric matrix

  for (I = N - 1; I >= 1; I--) {
    // generate random reflection

    dlarnv(3, ISEED, N - I + 1, WORK);
    WN = dnrm2(N - I + 1, WORK, 1);
    WA = sign(WN, WORK[1]).toDouble();
    if (WN == ZERO) {
      TAU = ZERO;
    } else {
      WB = WORK[1] + WA;
      dscal(N - I, ONE / WB, WORK(2), 1);
      WORK[1] = ONE;
      TAU = WB / WA;
    }

    // apply random reflection to A(i:n,i:n) from the left
    // and the right

    // compute  y := tau * A * u

    dsymv('Lower', N - I + 1, TAU, A(I, I), LDA, WORK, 1, ZERO, WORK(N + 1), 1);

    // compute  v := y - 1/2 * tau * ( y, u ) * u

    ALPHA = -HALF * TAU * ddot(N - I + 1, WORK(N + 1), 1, WORK, 1);
    daxpy(N - I + 1, ALPHA, WORK, 1, WORK(N + 1), 1);

    // apply the transformation as a rank-2 update to A(i:n,i:n)

    dsyr2('Lower', N - I + 1, -ONE, WORK, 1, WORK(N + 1), 1, A(I, I), LDA);
  }

  // Reduce number of subdiagonals to K

  for (I = 1; I <= N - 1 - K; I++) {
    // generate reflection to annihilate A(k+i+1:n,i)

    WN = dnrm2(N - K - I + 1, A(K + I, I).asArray(), 1);
    WA = sign(WN, A[K + I][I]).toDouble();
    if (WN == ZERO) {
      TAU = ZERO;
    } else {
      WB = A[K + I][I] + WA;
      dscal(N - K - I, ONE / WB, A(K + I + 1, I).asArray(), 1);
      A[K + I][I] = ONE;
      TAU = WB / WA;
    }

    // apply reflection to A(k+i:n,i+1:k+i-1) from the left

    dgemv('Transpose', N - K - I + 1, K - 1, ONE, A(K + I, I + 1), LDA,
        A(K + I, I).asArray(), 1, ZERO, WORK, 1);
    dger(N - K - I + 1, K - 1, -TAU, A(K + I, I).asArray(), 1, WORK, 1,
        A(K + I, I + 1), LDA);

    // apply reflection to A(k+i:n,k+i:n) from the left and the right

    // compute  y := tau * A * u

    dsymv('Lower', N - K - I + 1, TAU, A(K + I, K + I), LDA,
        A(K + I, I).asArray(), 1, ZERO, WORK, 1);

    // compute  v := y - 1/2 * tau * ( y, u ) * u

    ALPHA =
        -HALF * TAU * ddot(N - K - I + 1, WORK, 1, A(K + I, I).asArray(), 1);
    daxpy(N - K - I + 1, ALPHA, A(K + I, I).asArray(), 1, WORK, 1);

    // apply symmetric rank-2 update to A(k+i:n,k+i:n)

    dsyr2('Lower', N - K - I + 1, -ONE, A(K + I, I).asArray(), 1, WORK, 1,
        A(K + I, K + I), LDA);

    A[K + I][I] = -WA;
    for (J = K + I + 1; J <= N; J++) {
      A[J][I] = ZERO;
    }
  }

  // Store full symmetric matrix

  for (J = 1; J <= N; J++) {
    for (I = J + 1; I <= N; I++) {
      A[J][I] = A[I][J];
    }
  }
}
