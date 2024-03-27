import 'dart:math';

import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/zherk.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zgeqrf.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlansy.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zungqr.dart';

import 'common.dart';

void zqrt01(
  final int M,
  final int N,
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
  const ROGUE = Complex(-1.0e+10, -1.0e+10);
  final INFO = Box(0);

  final MINMN = min(M, N);
  final EPS = dlamch('Epsilon');

  // Copy the matrix A to the array AF.

  zlacpy('Full', M, N, A, LDA, AF, LDA);

  // Factorize the matrix A in the array AF.

  srnamc.SRNAMT = 'ZGEQRF';
  zgeqrf(M, N, AF, LDA, TAU, WORK, LWORK, INFO);

  // Copy details of Q

  zlaset('Full', M, M, ROGUE, ROGUE, Q, LDA);
  zlacpy('Lower', M - 1, N, AF(2, 1), LDA, Q(2, 1), LDA);

  // Generate the m-by-m matrix Q

  srnamc.SRNAMT = 'ZUNGQR';
  zungqr(M, M, MINMN, Q, LDA, TAU, WORK, LWORK, INFO);

  // Copy R

  zlaset('Full', M, N, Complex.zero, Complex.zero, R, LDA);
  zlacpy('Upper', M, N, AF, LDA, R, LDA);

  // Compute R - Q'*A

  zgemm('Conjugate transpose', 'No transpose', M, N, M, -Complex.one, Q, LDA, A,
      LDA, Complex.one, R, LDA);

  // Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .

  final ANORM = zlange('1', M, N, A, LDA, RWORK);
  var RESID = zlange('1', M, N, R, LDA, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / max(1, M)) / ANORM) / EPS;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute I - Q'*Q

  zlaset('Full', M, M, Complex.zero, Complex.one, R, LDA);
  zherk('Upper', 'Conjugate transpose', M, M, -ONE, Q, LDA, ONE, R, LDA);

  // Compute norm( I - Q'*Q ) / ( M * EPS ) .

  RESID = zlansy('1', 'Upper', M, R, LDA, RWORK);

  RESULT[2] = (RESID / max(1, M)) / EPS;
}
