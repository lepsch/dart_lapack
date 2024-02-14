import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlaswp.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dgesc2(
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> RHS_,
  final Array<int> IPIV_,
  final Array<int> JPIV_,
  final Box<double> SCALE,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final RHS = RHS_.dim();
  final IPIV = IPIV_.dim();
  final JPIV = JPIV_.dim();
  const ONE = 1.0, TWO = 2.0;
  int I, J;
  double
      // BIGNUM,
      EPS,
      SMLNUM,
      TEMP;

  // Set constant to control overflow

  EPS = dlamch('P');
  SMLNUM = dlamch('S') / EPS;
  // BIGNUM = ONE / SMLNUM;

  // Apply permutations IPIV to RHS

  dlaswp(1, RHS.asMatrix(LDA), LDA, 1, N - 1, IPIV, 1);

  // Solve for L part

  for (I = 1; I <= N - 1; I++) {
    for (J = I + 1; J <= N; J++) {
      RHS[J] = RHS[J] - A[J][I] * RHS[I];
    }
  }

  // Solve for U part

  SCALE.value = ONE;

  // Check for scaling

  I = idamax(N, RHS, 1);
  if (TWO * SMLNUM * (RHS[I]).abs() > (A[N][N]).abs()) {
    TEMP = (ONE / TWO) / (RHS[I]).abs();
    dscal(N, TEMP, RHS(1), 1);
    SCALE.value = SCALE.value * TEMP;
  }

  for (I = N; I >= 1; I--) {
    TEMP = ONE / A[I][I];
    RHS[I] = RHS[I] * TEMP;
    for (J = I + 1; J <= N; J++) {
      RHS[I] = RHS[I] - RHS[J] * (A[I][J] * TEMP);
    }
  }

  // Apply permutations JPIV to the solution (RHS)

  dlaswp(1, RHS.asMatrix(LDA), LDA, 1, N - 1, JPIV, -1);
}
