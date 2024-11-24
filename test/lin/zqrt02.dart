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
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlansy.dart';
import 'package:dart_lapack/src/zlaset.dart';
import 'package:dart_lapack/src/zungqr.dart';

import 'common.dart';

void zqrt02(
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
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDA);
  final Q = Q_.having(ld: LDA);
  final R = R_.having(ld: LDA);
  final TAU = TAU_.having();
  final RESULT = RESULT_.having();
  final WORK = WORK_.having(length: LWORK);
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  const ROGUE = Complex(-1.0e+10, -1.0e+10);
  final INFO = Box(0);

  final EPS = dlamch('Epsilon');

  // Copy the first k columns of the factorization to the array Q

  zlaset('Full', M, N, ROGUE, ROGUE, Q, LDA);
  zlacpy('Lower', M - 1, K, AF(2, 1), LDA, Q(2, 1), LDA);

  // Generate the first n columns of the matrix Q

  srnamc.SRNAMT = 'ZUNGQR';
  zungqr(M, N, K, Q, LDA, TAU, WORK, LWORK, INFO);

  // Copy R(1:n,1:k)

  zlaset('Full', N, K, Complex.zero, Complex.zero, R, LDA);
  zlacpy('Upper', N, K, AF, LDA, R, LDA);

  // Compute R(1:n,1:k) - Q(1:m,1:n)' * A(1:m,1:k)

  zgemm('Conjugate transpose', 'No transpose', N, K, M, -Complex.one, Q, LDA, A,
      LDA, Complex.one, R, LDA);

  // Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .

  final ANORM = zlange('1', M, K, A, LDA, RWORK);
  var RESID = zlange('1', N, K, R, LDA, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / max(1, M)) / ANORM) / EPS;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute I - Q'*Q

  zlaset('Full', N, N, Complex.zero, Complex.one, R, LDA);
  zherk('Upper', 'Conjugate transpose', N, M, -ONE, Q, LDA, ONE, R, LDA);

  // Compute norm( I - Q'*Q ) / ( M * EPS ) .

  RESID = zlansy('1', 'Upper', N, R, LDA, RWORK);

  RESULT[2] = (RESID / max(1, M)) / EPS;
}
