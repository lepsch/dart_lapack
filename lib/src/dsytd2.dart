import 'dart:math';

import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dsymv.dart';
import 'package:lapack/src/blas/dsyr2.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlarfg.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dsytd2(
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> D_,
  final Array<double> E_,
  final Array<double> TAU_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final D = D_.having();
  final E = E_.having();
  final TAU = TAU_.having();
  const ONE = 1.0, ZERO = 0.0, HALF = 1.0 / 2.0;
  bool UPPER;
  int I;
  double ALPHA;
  final TAUI = Box(0.0);

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
    xerbla('DSYTD2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N <= 0) return;

  if (UPPER) {
    // Reduce the upper triangle of A

    for (I = N - 1; I >= 1; I--) {
      // Generate elementary reflector H(i) = I - tau * v * v**T
      // to annihilate A(1:i-1,i+1)

      dlarfg(I, A.box(I, I + 1), A(1, I + 1).asArray(), 1, TAUI);
      E[I] = A[I][I + 1];

      if (TAUI.value != ZERO) {
        // Apply H(i) from both sides to A(1:i,1:i)

        A[I][I + 1] = ONE;

        // Compute  x := tau * A * v  storing x in TAU(1:i)

        dsymv(UPLO, I, TAUI.value, A, LDA, A(1, I + 1).asArray(), 1, ZERO, TAU,
            1);

        // Compute  w := x - 1/2 * tau * (x**T * v) * v

        ALPHA = -HALF * TAUI.value * ddot(I, TAU, 1, A(1, I + 1).asArray(), 1);
        daxpy(I, ALPHA, A(1, I + 1).asArray(), 1, TAU, 1);

        // Apply the transformation as a rank-2 update:
        //    A := A - v * w**T - w * v**T

        dsyr2(UPLO, I, -ONE, A(1, I + 1).asArray(), 1, TAU, 1, A, LDA);

        A[I][I + 1] = E[I];
      }
      D[I + 1] = A[I + 1][I + 1];
      TAU[I] = TAUI.value;
    }
    D[1] = A[1][1];
  } else {
    // Reduce the lower triangle of A

    for (I = 1; I <= N - 1; I++) {
      // Generate elementary reflector H(i) = I - tau * v * v**T
      // to annihilate A(i+2:n,i)

      dlarfg(N - I, A.box(I + 1, I), A(min(I + 2, N), I).asArray(), 1, TAUI);
      E[I] = A[I + 1][I];

      if (TAUI.value != ZERO) {
        // Apply H(i) from both sides to A(i+1:n,i+1:n)

        A[I + 1][I] = ONE;

        // Compute  x := tau * A * v  storing y in TAU(i:n-1)

        dsymv(UPLO, N - I, TAUI.value, A(I + 1, I + 1), LDA,
            A(I + 1, I).asArray(), 1, ZERO, TAU(I), 1);

        // Compute  w := x - 1/2 * tau * (x**T * v) * v

        ALPHA = -HALF *
            TAUI.value *
            ddot(N - I, TAU(I), 1, A(I + 1, I).asArray(), 1);
        daxpy(N - I, ALPHA, A(I + 1, I).asArray(), 1, TAU(I), 1);

        // Apply the transformation as a rank-2 update:
        //    A := A - v * w**T - w * v**T

        dsyr2(UPLO, N - I, -ONE, A(I + 1, I).asArray(), 1, TAU(I), 1,
            A(I + 1, I + 1), LDA);

        A[I + 1][I] = E[I];
      }
      D[I] = A[I][I];
      TAU[I] = TAUI.value;
    }
    D[N] = A[N][N];
  }
}
