import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dptcon(
  final int N,
  final Array<double> D_,
  final Array<double> E_,
  final double ANORM,
  final Box<double> RCOND,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final E = E_.having();
  final WORK = WORK_.having();

  const ONE = 1.0, ZERO = 0.0;
  int I, IX;
  double AINVNM;

  // Test the input arguments.

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
  } else if (ANORM < ZERO) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('DPTCON', -INFO.value);
    return;
  }

  // Quick return if possible

  RCOND.value = ZERO;
  if (N == 0) {
    RCOND.value = ONE;
    return;
  } else if (ANORM == ZERO) {
    return;
  }

  // Check that D(1:N) is positive.

  for (I = 1; I <= N; I++) {
    if (D[I] <= ZERO) return;
  }

  // Solve M(A) * x = e, where M(A) = (m(i,j)) is given by

  // m(i,j) =  abs(A(i,j)), i = j,
  // m(i,j) = -abs(A(i,j)), i != j,

  // and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**T.

  // Solve M(L) * x = e.

  WORK[1] = ONE;
  for (I = 2; I <= N; I++) {
    WORK[I] = ONE + WORK[I - 1] * E[I - 1].abs();
  }

  // Solve D * M(L)**T * x = b.

  WORK[N] /= D[N];
  for (I = N - 1; I >= 1; I--) {
    WORK[I] = WORK[I] / D[I] + WORK[I + 1] * E[I].abs();
  }

  // Compute AINVNM = max(x(i)), 1<=i<=n.

  IX = idamax(N, WORK, 1);
  AINVNM = WORK[IX].abs();

  // Compute the reciprocal condition number.

  if (AINVNM != ZERO) RCOND.value = (ONE / AINVNM) / ANORM;
}
