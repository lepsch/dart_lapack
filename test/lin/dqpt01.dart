import 'dart:math';

import 'package:lapack/lapack.dart';

double dqpt01(
  final int M,
  final int N,
  final int K,
  final Matrix<double> A_,
  final Matrix<double> AF_,
  final int LDA,
  final Array<double> TAU_,
  final Array<int> JPVT_,
  final Array<double> WORK_,
  final int LWORK,
) {
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDA);
  final TAU = TAU_.having();
  final JPVT = JPVT_.having();
  final WORK = WORK_.having(length: LWORK);

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  final RWORK = Array<double>(1);
  final INFO = Box(0);

  // Test if there is enough workspace

  if (LWORK < M * N + N) {
    xerbla('DQPT01', 10);
    return ZERO;
  }

  // Quick return if possible

  if (M <= 0 || N <= 0) return ZERO;

  final NORMA = dlange('One-norm', M, N, A, LDA, RWORK);

  for (var J = 1; J <= K; J++) {
    // Copy the upper triangular part of the factor R stored
    // in AF(1:K,1:K) into the work array WORK.

    for (var I = 1; I <= min(J, M); I++) {
      WORK[(J - 1) * M + I] = AF[I][J];
    }

    // Zero out the elements below the diagonal in the work array.

    for (var I = J + 1; I <= M; I++) {
      WORK[(J - 1) * M + I] = ZERO;
    }
  }

  // Copy columns (K+1,N) from AF into the work array WORK.
  // AF(1:K,K+1:N) contains the rectangular block of the upper trapezoidal
  // factor R, AF(K+1:M,K+1:N) contains the partially updated residual
  // matrix of R.

  for (var J = K + 1; J <= N; J++) {
    dcopy(M, AF(1, J).asArray(), 1, WORK((J - 1) * M + 1), 1);
  }

  dormqr('Left', 'No transpose', M, N, K, AF, LDA, TAU, WORK.asMatrix(), M,
      WORK(M * N + 1), LWORK - M * N, INFO);

  for (var J = 1; J <= N; J++) {
    // Compare J-th column of QR and JPVT(J)-th column of A.

    daxpy(M, -ONE, A(1, JPVT[J]).asArray(), 1, WORK((J - 1) * M + 1), 1);
  }

  final result = dlange('One-norm', M, N, WORK.asMatrix(), M, RWORK) /
      (max(M, N) * dlamch('Epsilon'));
  return NORMA != ZERO ? result / NORMA : result;
}
