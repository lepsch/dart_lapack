import 'dart:math';

import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/zherk.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlansy.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zungql.dart';

import 'common.dart';

void zqlt02(
  final int M,
  final int N,
  final int K,
  final Matrix<Complex> A_,
  final Matrix<Complex> AF_,
  final Matrix<Complex> Q_,
  final Matrix<Complex> L_,
  final int LDA,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<double> RESULT_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDA);
  final Q = Q_.having(ld: LDA);
  final L = L_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having(length: LWORK);
  final RWORK = RWORK_.having();
  final RESULT = RESULT_.having();
  const ZERO = 0.0, ONE = 1.0;
  const ROGUE = Complex(-1.0e+10, -1.0e+10);
  final INFO = Box(0);

  // Quick return if possible

  if (M == 0 || N == 0 || K == 0) {
    RESULT[1] = ZERO;
    RESULT[2] = ZERO;
    return;
  }

  final EPS = dlamch('Epsilon');

  // Copy the last k columns of the factorization to the array Q

  zlaset('Full', M, N, ROGUE, ROGUE, Q, LDA);
  if (K < M) {
    zlacpy('Full', M - K, K, AF(1, N - K + 1), LDA, Q(1, N - K + 1), LDA);
  }
  if (K > 1) {
    zlacpy('Upper', K - 1, K - 1, AF(M - K + 1, N - K + 2), LDA,
        Q(M - K + 1, N - K + 2), LDA);
  }

  // Generate the last n columns of the matrix Q

  srnamc.SRNAMT = 'ZUNGQL';
  zungql(M, N, K, Q, LDA, TAU(N - K + 1), WORK, LWORK, INFO);

  // Copy L(m-n+1:m,n-k+1:n)

  zlaset(
      'Full', N, K, Complex.zero, Complex.zero, L(M - N + 1, N - K + 1), LDA);
  zlacpy('Lower', K, K, AF(M - K + 1, N - K + 1), LDA, L(M - K + 1, N - K + 1),
      LDA);

  // Compute L(m-n+1:m,n-k+1:n) - Q(1:m,m-n+1:m)' * A(1:m,n-k+1:n)

  zgemm('Conjugate transpose', 'No transpose', N, K, M, -Complex.one, Q, LDA,
      A(1, N - K + 1), LDA, Complex.one, L(M - N + 1, N - K + 1), LDA);

  // Compute norm( L - Q'*A ) / ( M * norm(A) * EPS ) .

  final ANORM = zlange('1', M, K, A(1, N - K + 1), LDA, RWORK);
  var RESID = zlange('1', N, K, L(M - N + 1, N - K + 1), LDA, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / (max(1, M)).toDouble()) / ANORM) / EPS;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute I - Q'*Q

  zlaset('Full', N, N, Complex.zero, Complex.one, L, LDA);
  zherk('Upper', 'Conjugate transpose', N, M, -ONE, Q, LDA, ONE, L, LDA);

  // Compute norm( I - Q'*Q ) / ( M * EPS ) .

  RESID = zlansy('1', 'Upper', N, L, LDA, RWORK);

  RESULT[2] = (RESID / (max(1, M)).toDouble()) / EPS;
}
