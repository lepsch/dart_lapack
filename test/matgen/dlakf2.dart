import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/matrix.dart';

void dlakf2(
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final Matrix<double> D_,
  final Matrix<double> E_,
  final Matrix<double> Z_,
  final int LDZ,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final B = B_.dim(LDA);
  final D = D_.dim(LDA);
  final E = E_.dim(LDA);
  final Z = Z_.dim(LDZ);
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
