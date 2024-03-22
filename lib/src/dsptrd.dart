import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dspmv.dart';
import 'package:lapack/src/blas/dspr2.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlarfg.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dsptrd(
  final String UPLO,
  final int N,
  final Array<double> AP_,
  final Array<double> D_,
  final Array<double> E_,
  final Array<double> TAU_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.having();
  final D = D_.having();
  final E = E_.having();
  final TAU = TAU_.having();
  const ONE = 1.0, ZERO = 0.0, HALF = 1.0 / 2.0;
  bool UPPER;
  int I, I1, I1I1, II;
  double ALPHA;
  final TAUI = Box(0.0);

  // Test the input parameters

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  }
  if (INFO.value != 0) {
    xerbla('DSPTRD', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N <= 0) return;

  if (UPPER) {
    // Reduce the upper triangle of A.
    // I1 is the index in AP of A[1,I+1].

    I1 = N * (N - 1) ~/ 2 + 1;
    for (I = N - 1; I >= 1; I--) {
      // Generate elementary reflector H(i) = I - tau * v * v**T
      // to annihilate A[1:i-1,i+1]

      dlarfg(I, AP(I1 + I - 1), AP(I1), 1, TAUI);
      E[I] = AP[I1 + I - 1];

      if (TAUI.value != ZERO) {
        // Apply H(i) from both sides to A[1:i,1:i]

        AP[I1 + I - 1] = ONE;

        // Compute  y := tau * A * v  storing y in TAU[1:i]

        dspmv(UPLO, I, TAUI.value, AP, AP(I1), 1, ZERO, TAU, 1);

        // Compute  w := y - 1/2 * tau * (y**T *v) * v

        ALPHA = -HALF * TAUI.value * ddot(I, TAU, 1, AP(I1), 1);
        daxpy(I, ALPHA, AP(I1), 1, TAU, 1);

        // Apply the transformation as a rank-2 update:
        // A := A - v * w**T - w * v**T

        dspr2(UPLO, I, -ONE, AP(I1), 1, TAU, 1, AP);

        AP[I1 + I - 1] = E[I];
      }
      D[I + 1] = AP[I1 + I];
      TAU[I] = TAUI.value;
      I1 -= I;
    }
    D[1] = AP[1];
  } else {
    // Reduce the lower triangle of A. II is the index in AP of
    // A[i,i] and I1I1 is the index of A[i+1,i+1].

    II = 1;
    for (I = 1; I <= N - 1; I++) {
      I1I1 = II + N - I + 1;

      // Generate elementary reflector H(i) = I - tau * v * v**T
      // to annihilate A[i+2:n,i]

      dlarfg(N - I, AP(II + 1), AP(II + 2), 1, TAUI);
      E[I] = AP[II + 1];

      if (TAUI.value != ZERO) {
        // Apply H(i) from both sides to A[i+1:n,i+1:n]

        AP[II + 1] = ONE;

        // Compute  y := tau * A * v  storing y in TAU[i:n-1]

        dspmv(
            UPLO, N - I, TAUI.value, AP(I1I1), AP(II + 1), 1, ZERO, TAU(I), 1);

        // Compute  w := y - 1/2 * tau * (y**T *v) * v

        ALPHA = -HALF * TAUI.value * ddot(N - I, TAU(I), 1, AP(II + 1), 1);
        daxpy(N - I, ALPHA, AP(II + 1), 1, TAU(I), 1);

        // Apply the transformation as a rank-2 update:
        // A := A - v * w**T - w * v**T

        dspr2(UPLO, N - I, -ONE, AP(II + 1), 1, TAU(I), 1, AP(I1I1));

        AP[II + 1] = E[I];
      }
      D[I] = AP[II];
      TAU[I] = TAUI.value;
      II = I1I1;
    }
    D[N] = AP[II];
  }
}
