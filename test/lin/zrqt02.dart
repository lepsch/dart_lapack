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
import 'package:lapack/src/zungrq.dart';

import 'common.dart';

void zrqt02(
  final int M,
  final int N,
  final int K,
  final Matrix<Complex> A_,
  final Matrix<Complex> AF_,
  final Matrix<Complex> Q_,
  final Matrix<Complex> R_,
  final int LDA,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<double> RESULT_,
) {
  final A = A_.having();
  final AF = AF_.having();
  final Q = Q_.having();
  final R = R_.having();
  final TAU = TAU_.having();
  final WORK = WORK_.having(length: LWORK);
  final RWORK = RWORK_.having();
  final RESULT = RESULT_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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

  // Copy the last k rows of the factorization to the array Q

  zlaset('Full', M, N, ROGUE, ROGUE, Q, LDA);
  if (K < N) {
    zlacpy('Full', K, N - K, AF(M - K + 1, 1), LDA, Q(M - K + 1, 1), LDA);
  }
  if (K > 1) {
    zlacpy('Lower', K - 1, K - 1, AF(M - K + 2, N - K + 1), LDA,
        Q(M - K + 2, N - K + 1), LDA);
  }

  // Generate the last n rows of the matrix Q

  srnamc.SRNAMT = 'ZUNGRQ';
  zungrq(M, N, K, Q, LDA, TAU(M - K + 1), WORK, LWORK, INFO);

  // Copy R(m-k+1:m,n-m+1:n)

  zlaset(
      'Full', K, M, Complex.zero, Complex.zero, R(M - K + 1, N - M + 1), LDA);
  zlacpy('Upper', K, K, AF(M - K + 1, N - K + 1), LDA, R(M - K + 1, N - K + 1),
      LDA);

  // Compute R(m-k+1:m,n-m+1:n) - A(m-k+1:m,1:n) * Q(n-m+1:n,1:n)'

  zgemm('No transpose', 'Conjugate transpose', K, M, N, Complex(-ONE),
      A(M - K + 1, 1), LDA, Q, LDA, Complex.one, R(M - K + 1, N - M + 1), LDA);

  // Compute norm( R - A*Q' ) / ( N * norm(A) * EPS ) .

  final ANORM = zlange('1', K, N, A(M - K + 1, 1), LDA, RWORK);
  var RESID = zlange('1', K, M, R(M - K + 1, N - M + 1), LDA, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / max(1, N)) / ANORM) / EPS;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute I - Q*Q'

  zlaset('Full', M, M, Complex.zero, Complex.one, R, LDA);
  zherk('Upper', 'No transpose', M, N, -ONE, Q, LDA, ONE, R, LDA);

  // Compute norm( I - Q*Q' ) / ( N * EPS ) .

  RESID = zlansy('1', 'Upper', M, R, LDA, RWORK);

  RESULT[2] = (RESID / max(1, N)) / EPS;
}
