import 'dart:math';

import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dger.dart';
import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlarnv.dart';
import 'package:lapack/src/f2c/sign.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlarge(
  final int N,
  final Matrix<double> A,
  final int LDA,
  final Array<int> ISEED,
  final Array<double> WORK,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  int I;
  double TAU, WA, WB, WN;

  // Test the input arguments

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
  } else if (LDA < max(1, N)) {
    INFO.value = -3;
  }
  if (INFO.value < 0) {
    xerbla('DLARGE', -INFO.value);
    return;
  }

  // pre- and post-multiply A by random orthogonal matrix

  for (I = N; I >= 1; I--) {
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

    // multiply A[i:n][1:n] by random reflection from the left

    dgemv('Transpose', N - I + 1, N, ONE, A(I, 1), LDA, WORK, 1, ZERO,
        WORK(N + 1), 1);
    dger(N - I + 1, N, -TAU, WORK, 1, WORK(N + 1), 1, A(I, 1), LDA);

    // multiply A[1:n][i:n] by random reflection from the right

    dgemv('No transpose', N, N - I + 1, ONE, A(1, I), LDA, WORK, 1, ZERO,
        WORK(N + 1), 1);
    dger(N, N - I + 1, -TAU, WORK(N + 1), 1, WORK, 1, A(1, I), LDA);
  }
}
