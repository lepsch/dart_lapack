import 'package:lapack/lapack.dart';

void zgeqls(
  final int M,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> TAU_,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final TAU = TAU_.having();
  final WORK = WORK_.having(length: LWORK);

  // Test the input arguments.

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0 || N > M) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, M)) {
    INFO.value = -5;
  } else if (LDB < max(1, M)) {
    INFO.value = -8;
  } else if (LWORK < 1 || LWORK < NRHS && M > 0 && N > 0) {
    INFO.value = -10;
  }
  if (INFO.value != 0) {
    xerbla('ZGEQLS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0 || M == 0) return;

  // B := Q' * B

  zunmql('Left', 'Conjugate transpose', M, NRHS, N, A, LDA, TAU, B, LDB, WORK,
      LWORK, INFO);

  // Solve L*X = B(m-n+1:m,:)

  ztrsm('Left', 'Lower', 'No transpose', 'Non-unit', N, NRHS, Complex.one,
      A(M - N + 1, 1), LDA, B(M - N + 1, 1), LDB);
}
