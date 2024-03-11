import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dsyrk.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgeqrfp.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dorgqr.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

import 'common.dart';

void dqrt01p(
  final int M,
  final int N,
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
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDA);
  final Q = Q_.having(ld: LDA);
  final R = R_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having(length: LWORK);
  final RWORK = RWORK_.having();
  final RESULT = RESULT_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  const ROGUE = -1.0e+10;
  final INFO = Box(0);

  final MINMN = min(M, N);
  final EPS = dlamch('Epsilon');

  // Copy the matrix A to the array AF.

  dlacpy('Full', M, N, A, LDA, AF, LDA);

  // Factorize the matrix A in the array AF.

  srnamc.SRNAMT = 'DGEQRFP';
  dgeqrfp(M, N, AF, LDA, TAU, WORK, LWORK, INFO);

  // Copy details of Q

  dlaset('Full', M, M, ROGUE, ROGUE, Q, LDA);
  dlacpy('Lower', M - 1, N, AF(2, 1), LDA, Q(2, 1), LDA);

  // Generate the m-by-m matrix Q

  srnamc.SRNAMT = 'DORGQR';
  dorgqr(M, M, MINMN, Q, LDA, TAU, WORK, LWORK, INFO);

  // Copy R

  dlaset('Full', M, N, ZERO, ZERO, R, LDA);
  dlacpy('Upper', M, N, AF, LDA, R, LDA);

  // Compute R - Q'*A

  dgemm(
      'Transpose', 'No transpose', M, N, M, -ONE, Q, LDA, A, LDA, ONE, R, LDA);

  // Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .

  final ANORM = dlange('1', M, N, A, LDA, RWORK);
  var RESID = dlange('1', M, N, R, LDA, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / (max(1, M)).toDouble()) / ANORM) / EPS;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute I - Q'*Q

  dlaset('Full', M, M, ZERO, ONE, R, LDA);
  dsyrk('Upper', 'Transpose', M, M, -ONE, Q, LDA, ONE, R, LDA);

  // Compute norm( I - Q'*Q ) / ( M * EPS ) .

  RESID = dlansy('1', 'Upper', M, R, LDA, RWORK);

  RESULT[2] = (RESID / (max(1, M)).toDouble()) / EPS;
}
