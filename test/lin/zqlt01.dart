import 'dart:math';

import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/zherk.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zgeqlf.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlansy.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zungql.dart';

import 'common.dart';

void zqlt01(
  final int M,
  final int N,
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

  final MINMN = min(M, N);
  final EPS = dlamch('Epsilon');

  // Copy the matrix A to the array AF.

  zlacpy('Full', M, N, A, LDA, AF, LDA);

  // Factorize the matrix A in the array AF.

  srnamc.SRNAMT = 'ZGEQLF';
  zgeqlf(M, N, AF, LDA, TAU, WORK, LWORK, INFO);

  // Copy details of Q

  zlaset('Full', M, M, ROGUE, ROGUE, Q, LDA);
  if (M >= N) {
    if (N < M && N > 0) zlacpy('Full', M - N, N, AF, LDA, Q(1, M - N + 1), LDA);
    if (N > 1) {
      zlacpy('Upper', N - 1, N - 1, AF(M - N + 1, 2), LDA,
          Q(M - N + 1, M - N + 2), LDA);
    }
  } else {
    if (M > 1) {
      zlacpy('Upper', M - 1, M - 1, AF(1, N - M + 2), LDA, Q(1, 2), LDA);
    }
  }

  // Generate the m-by-m matrix Q

  srnamc.SRNAMT = 'ZUNGQL';
  zungql(M, M, MINMN, Q, LDA, TAU, WORK, LWORK, INFO);

  // Copy L

  zlaset('Full', M, N, Complex.zero, Complex.zero, L, LDA);
  if (M >= N) {
    if (N > 0) {
      zlacpy('Lower', N, N, AF(M - N + 1, 1), LDA, L(M - N + 1, 1), LDA);
    }
  } else {
    if (N > M && M > 0) zlacpy('Full', M, N - M, AF, LDA, L, LDA);
    if (M > 0) {
      zlacpy('Lower', M, M, AF(1, N - M + 1), LDA, L(1, N - M + 1), LDA);
    }
  }

  // Compute L - Q'*A

  zgemm('Conjugate transpose', 'No transpose', M, N, M, -Complex.one, Q, LDA, A,
      LDA, Complex.one, L, LDA);

  // Compute norm( L - Q'*A ) / ( M * norm(A) * EPS ) .

  final ANORM = zlange('1', M, N, A, LDA, RWORK);
  var RESID = zlange('1', M, N, L, LDA, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / max(1, M)) / ANORM) / EPS;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute I - Q'*Q

  zlaset('Full', M, M, Complex.zero, Complex.one, L, LDA);
  zherk('Upper', 'Conjugate transpose', M, M, -ONE, Q, LDA, ONE, L, LDA);

  // Compute norm( I - Q'*Q ) / ( M * EPS ) .

  RESID = zlansy('1', 'Upper', M, L, LDA, RWORK);

  RESULT[2] = (RESID / max(1, M)) / EPS;
}
