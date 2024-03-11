import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dsyrk.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dorgrq.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

import 'common.dart';

void drqt02(
  final int M,
  final int N,
  final int K,
  final Matrix<double> A_,
  final Matrix<double> AF_,
  final Matrix<double> Q_,
  final Matrix<double> R_,
  final int LDA,
  final Array<double> TAU_,
  final Array<double> WORK_,
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
  final R = R_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having(length: LWORK);
  final RWORK = RWORK_.having();
  final RESULT = RESULT_.having();
  const ZERO = 0.0, ONE = 1.0;
  const ROGUE = -1.0e+10;
  final INFO = Box(0);

  // Quick return if possible

  if (M == 0 || N == 0 || K == 0) {
    RESULT[1] = ZERO;
    RESULT[2] = ZERO;
    return;
  }

  final EPS = dlamch('Epsilon');

  // Copy the last k rows of the factorization to the array Q

  dlaset('Full', M, N, ROGUE, ROGUE, Q, LDA);
  if (K < N) {
    dlacpy('Full', K, N - K, AF(M - K + 1, 1), LDA, Q(M - K + 1, 1), LDA);
  }
  if (K > 1) {
    dlacpy('Lower', K - 1, K - 1, AF(M - K + 2, N - K + 1), LDA,
        Q(M - K + 2, N - K + 1), LDA);
  }

  // Generate the last n rows of the matrix Q

  srnamc.SRNAMT = 'DORGRQ';
  dorgrq(M, N, K, Q, LDA, TAU(M - K + 1), WORK, LWORK, INFO);

  // Copy R(m-k+1:m,n-m+1:n)

  dlaset('Full', K, M, ZERO, ZERO, R(M - K + 1, N - M + 1), LDA);
  dlacpy('Upper', K, K, AF(M - K + 1, N - K + 1), LDA, R(M - K + 1, N - K + 1),
      LDA);

  // Compute R(m-k+1:m,n-m+1:n) - A(m-k+1:m,1:n) * Q(n-m+1:n,1:n)'

  dgemm('No transpose', 'Transpose', K, M, N, -ONE, A(M - K + 1, 1), LDA, Q,
      LDA, ONE, R(M - K + 1, N - M + 1), LDA);

  // Compute norm( R - A*Q' ) / ( N * norm(A) * EPS ) .

  final ANORM = dlange('1', K, N, A(M - K + 1, 1), LDA, RWORK);
  var RESID = dlange('1', K, M, R(M - K + 1, N - M + 1), LDA, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / max(1, N)) / ANORM) / EPS;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute I - Q*Q'

  dlaset('Full', M, M, ZERO, ONE, R, LDA);
  dsyrk('Upper', 'No transpose', M, N, -ONE, Q, LDA, ONE, R, LDA);

  // Compute norm( I - Q*Q' ) / ( N * EPS ) .

  RESID = dlansy('1', 'Upper', M, R, LDA, RWORK);

  RESULT[2] = (RESID / max(1, N)) / EPS;
}
