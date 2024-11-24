// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlarnv.dart';
import 'package:dart_lapack/src/zlaset.dart';
import 'package:dart_lapack/src/zungqr.dart';
import 'package:dart_lapack/src/zunmqr.dart';

import 'common.dart';

void zqrt03(
  final int M,
  final int N,
  final int K,
  final Matrix<Complex> AF_,
  final Matrix<Complex> C_,
  final Matrix<Complex> CC_,
  final Matrix<Complex> Q_,
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
  final AF = AF_.having(ld: LDA);
  final C = C_.having(ld: LDA);
  final CC = CC_.having(ld: LDA);
  final Q = Q_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having(length: LWORK);
  final RWORK = RWORK_.having();
  final RESULT = RESULT_.having();
  const ZERO = 0.0, ONE = 1.0;
  const ROGUE = Complex(-1.0e+10, -1.0e+10);
  final ISEED = Array.fromList([1988, 1989, 1990, 1991]);
  final INFO = Box(0);

  final EPS = dlamch('Epsilon');

  // Copy the first k columns of the factorization to the array Q

  zlaset('Full', M, M, ROGUE, ROGUE, Q, LDA);
  zlacpy('Lower', M - 1, K, AF(2, 1), LDA, Q(2, 1), LDA);

  // Generate the m-by-m matrix Q

  srnamc.SRNAMT = 'ZUNGQR';
  zungqr(M, M, K, Q, LDA, TAU, WORK, LWORK, INFO);

  for (var ISIDE = 1; ISIDE <= 2; ISIDE++) {
    final (SIDE, MC, NC) = ISIDE == 1 ? ('L', M, N) : ('R', N, M);

    // Generate MC by NC matrix C

    for (var J = 1; J <= NC; J++) {
      zlarnv(2, ISEED, MC, C(1, J).asArray());
    }
    var CNORM = zlange('1', MC, NC, C, LDA, RWORK);
    if (CNORM == ZERO) CNORM = ONE;

    for (var ITRANS = 1; ITRANS <= 2; ITRANS++) {
      final TRANS = ITRANS == 1 ? 'N' : 'C';

      // Copy C

      zlacpy('Full', MC, NC, C, LDA, CC, LDA);

      // Apply Q or Q' to C

      srnamc.SRNAMT = 'ZUNMQR';
      zunmqr(SIDE, TRANS, MC, NC, K, AF, LDA, TAU, CC, LDA, WORK, LWORK, INFO);

      // Form explicit product and subtract

      if (lsame(SIDE, 'L')) {
        zgemm(TRANS, 'No transpose', MC, NC, MC, -Complex.one, Q, LDA, C, LDA,
            Complex.one, CC, LDA);
      } else {
        zgemm('No transpose', TRANS, MC, NC, NC, -Complex.one, C, LDA, Q, LDA,
            Complex.one, CC, LDA);
      }

      // Compute error in the difference

      final RESID = zlange('1', MC, NC, CC, LDA, RWORK);
      RESULT[(ISIDE - 1) * 2 + ITRANS] = RESID / (max(1, M) * CNORM * EPS);
    }
  }
}
