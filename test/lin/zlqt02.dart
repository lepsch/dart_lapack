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
import 'package:lapack/src/zunglq.dart';

import 'common.dart';

void zlqt02(
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

  final EPS = dlamch('Epsilon');

  // Copy the first k rows of the factorization to the array Q

  zlaset('Full', M, N, ROGUE, ROGUE, Q, LDA);
  zlacpy('Upper', K, N - 1, AF(1, 2), LDA, Q(1, 2), LDA);

  // Generate the first n columns of the matrix Q

  srnamc.SRNAMT = 'ZUNGLQ';
  zunglq(M, N, K, Q, LDA, TAU, WORK, LWORK, INFO);

  // Copy L(1:k,1:m)

  zlaset('Full', K, M, Complex.zero, Complex.zero, L, LDA);
  zlacpy('Lower', K, M, AF, LDA, L, LDA);

  // Compute L(1:k,1:m) - A(1:k,1:n) * Q(1:m,1:n)'

  zgemm('No transpose', 'Conjugate transpose', K, M, N, Complex(-ONE), A, LDA,
      Q, LDA, Complex.one, L, LDA);

  // Compute norm( L - A*Q' ) / ( N * norm(A) * EPS ) .

  final ANORM = zlange('1', K, N, A, LDA, RWORK);
  var RESID = zlange('1', K, M, L, LDA, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / max(1, N)) / ANORM) / EPS;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute I - Q*Q'

  zlaset('Full', M, M, Complex.zero, Complex.one, L, LDA);
  zherk('Upper', 'No transpose', M, N, -ONE, Q, LDA, ONE, L, LDA);

  // Compute norm( I - Q*Q' ) / ( N * EPS ) .

  RESID = zlansy('1', 'Upper', M, L, LDA, RWORK);

  RESULT[2] = (RESID / max(1, N)) / EPS;
}
