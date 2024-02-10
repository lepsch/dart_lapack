import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/matrix.dart';

void dlakf2(
  final int M,
  final int N,
  final Matrix<double> A,
  final int LDA,
  final Matrix<double> B,
  final Matrix<double> D,
  final Matrix<double> E,
  final Matrix<double> Z,
  final int LDZ,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;
  int I, IK, J, JK, L, MN, MN2;

  // Initialize Z

  MN = M * N;
  MN2 = 2 * MN;
  dlaset('Full', MN2, MN2, ZERO, ZERO, Z, LDZ);

  IK = 1;
  for (L = 1; L <= N; L++) {
    // form kron(In, A)

    for (I = 1; I <= M; I++) {
      for (J = 1; J <= M; J++) {
        Z[IK + I - 1][IK + J - 1] = A[I][J];
      }
    }

    // form kron(In, D)

    for (I = 1; I <= M; I++) {
      for (J = 1; J <= M; J++) {
        Z[IK + MN + I - 1][IK + J - 1] = D[I][J];
      }
    }

    IK = IK + M;
  }

  IK = 1;
  for (L = 1; L <= N; L++) {
    JK = MN + 1;

    for (J = 1; J <= N; J++) {
      // form -kron(B', Im)

      for (I = 1; I <= M; I++) {
        Z[IK + I - 1][JK + I - 1] = -B[J][L];
      }

      // form -kron(E', Im)

      for (I = 1; I <= M; I++) {
        Z[IK + MN + I - 1][JK + I - 1] = -E[J][L];
      }

      JK = JK + M;
    }

    IK = IK + M;
  }
}
