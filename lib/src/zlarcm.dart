import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zlarcm(
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> C_,
  final int LDC,
  final Array<double> RWORK_,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final B = B_.dim(LDB);
  final C = C_.dim(LDC);
  final RWORK = RWORK_.dim();
  const ONE = 1.0, ZERO = 0.0;
  int I, J, L;

  // Quick return if possible.

  if ((M == 0) || (N == 0)) return;

  for (J = 1; J <= N; J++) {
    // 20
    for (I = 1; I <= M; I++) {
      // 10
      RWORK[(J - 1) * M + I] = (B[I][J]).toDouble();
    } // 10
  } // 20

  L = M * N + 1;
  dgemm('N', 'N', M, N, M, ONE, A, LDA, RWORK.asMatrix(M), M, ZERO,
      RWORK(L).asMatrix(M), M);
  for (J = 1; J <= N; J++) {
    // 40
    for (I = 1; I <= M; I++) {
      // 30
      C[I][J] = RWORK[L + (J - 1) * M + I - 1].toComplex();
    } // 30
  } // 40

  for (J = 1; J <= N; J++) {
    // 60
    for (I = 1; I <= M; I++) {
      // 50
      RWORK[(J - 1) * M + I] = B[I][J].imaginary;
    } // 50
  } // 60
  dgemm('N', 'N', M, N, M, ONE, A, LDA, RWORK.asMatrix(M), M, ZERO,
      RWORK(L).asMatrix(M), M);
  for (J = 1; J <= N; J++) {
    // 80
    for (I = 1; I <= M; I++) {
      // 70
      C[I][J] = Complex(C[I][J].toDouble(), RWORK[L + (J - 1) * M + I - 1]);
    } // 70
  } // 80
}
