import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zlacrm(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<Complex> C_,
  final int LDC,
  final Array<double> RWORK_,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final C = C_.having(ld: LDC);
  final RWORK = RWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  int I, J, L;

  // Quick return if possible.

  if ((M == 0) || (N == 0)) return;

  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      RWORK[(J - 1) * M + I] = (A[I][J]).toDouble();
    }
  }

  L = M * N + 1;
  dgemm('N', 'N', M, N, N, ONE, RWORK.asMatrix(M), M, B, LDB, ZERO,
      RWORK(L).asMatrix(M), M);
  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      C[I][J] = RWORK[L + (J - 1) * M + I - 1].toComplex();
    }
  }

  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      RWORK[(J - 1) * M + I] = A[I][J].imaginary;
    }
  }
  dgemm('N', 'N', M, N, N, ONE, RWORK.asMatrix(M), M, B, LDB, ZERO,
      RWORK(L).asMatrix(M), M);
  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      C[I][J] = Complex((C[I][J]).toDouble(), RWORK[L + (J - 1) * M + I - 1]);
    }
  }
}
