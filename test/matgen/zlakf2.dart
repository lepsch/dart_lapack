import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlaset.dart';

void zlakf2(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final Matrix<Complex> D_,
  final Matrix<Complex> E_,
  final Matrix<Complex> Z_,
  final int LDZ,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDA);
  final D = D_.having(ld: LDA);
  final E = E_.having(ld: LDA);
  final Z = Z_.having(ld: LDA);
  int I, IK, J, JK, L, MN, MN2;
  // ..
  // .. External Subroutines ..
  // EXTERNAL ZLASET

  // Initialize Z

  MN = M * N;
  MN2 = 2 * MN;
  zlaset('Full', MN2, MN2, Complex.zero, Complex.zero, Z, LDZ);

  IK = 1;
  for (L = 1; L <= N; L++) {
    // 50

    // form kron(In, A)

    for (I = 1; I <= M; I++) {
      // 20
      for (J = 1; J <= M; J++) {
        // 10
        Z[IK + I - 1][IK + J - 1] = A[I][J];
      } // 10
    } // 20

    // form kron(In, D)

    for (I = 1; I <= M; I++) {
      // 40
      for (J = 1; J <= M; J++) {
        // 30
        Z[IK + MN + I - 1][IK + J - 1] = D[I][J];
      } // 30
    } // 40

    IK = IK + M;
  } // 50

  IK = 1;
  for (L = 1; L <= N; L++) {
    // 90
    JK = MN + 1;

    for (J = 1; J <= N; J++) {
      // 80

      // form -kron(B', Im)

      for (I = 1; I <= M; I++) {
        // 60
        Z[IK + I - 1][JK + I - 1] = -B[J][L];
      } // 60

      // form -kron(E', Im)

      for (I = 1; I <= M; I++) {
        // 70
        Z[IK + MN + I - 1][JK + I - 1] = -E[J][L];
      } // 70

      JK = JK + M;
    } // 80

    IK = IK + M;
  } // 90
}
