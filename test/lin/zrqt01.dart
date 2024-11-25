// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/blas/zherk.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zgerqf.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlansy.dart';
import 'package:dart_lapack/src/zlaset.dart';
import 'package:dart_lapack/src/zungrq.dart';

import 'common.dart';

void zrqt01(
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

  srnamc.SRNAMT = 'ZGERQF';
  zgerqf(M, N, AF, LDA, TAU, WORK, LWORK, INFO);

  // Copy details of Q

  zlaset('Full', N, N, ROGUE, ROGUE, Q, LDA);
  if (M <= N) {
    if (M > 0 && M < N) zlacpy('Full', M, N - M, AF, LDA, Q(N - M + 1, 1), LDA);
    if (M > 1) {
      zlacpy('Lower', M - 1, M - 1, AF(2, N - M + 1), LDA,
          Q(N - M + 2, N - M + 1), LDA);
    }
  } else {
    if (N > 1) {
      zlacpy('Lower', N - 1, N - 1, AF(M - N + 2, 1), LDA, Q(2, 1), LDA);
    }
  }

  // Generate the n-by-n matrix Q

  srnamc.SRNAMT = 'ZUNGRQ';
  zungrq(N, N, MINMN, Q, LDA, TAU, WORK, LWORK, INFO);

  // Copy R

  zlaset('Full', M, N, Complex.zero, Complex.zero, R, LDA);
  if (M <= N) {
    if (M > 0) {
      zlacpy('Upper', M, M, AF(1, N - M + 1), LDA, R(1, N - M + 1), LDA);
    }
  } else {
    if (M > N && N > 0) zlacpy('Full', M - N, N, AF, LDA, R, LDA);
    if (N > 0) {
      zlacpy('Upper', N, N, AF(M - N + 1, 1), LDA, R(M - N + 1, 1), LDA);
    }
  }

  // Compute R - A*Q'

  zgemm('No transpose', 'Conjugate transpose', M, N, N, Complex(-ONE), A, LDA,
      Q, LDA, Complex.one, R, LDA);

  // Compute norm( R - Q'*A ) / ( N * norm(A) * EPS ) .

  final ANORM = zlange('1', M, N, A, LDA, RWORK);
  var RESID = zlange('1', M, N, R, LDA, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / max(1, N)) / ANORM) / EPS;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute I - Q*Q'

  zlaset('Full', N, N, Complex.zero, Complex.one, R, LDA);
  zherk('Upper', 'No transpose', N, N, -ONE, Q, LDA, ONE, R, LDA);

  // Compute norm( I - Q*Q' ) / ( N * EPS ) .

  RESID = zlansy('1', 'Upper', N, R, LDA, RWORK);

  RESULT[2] = (RESID / max(1, N)) / EPS;
}
