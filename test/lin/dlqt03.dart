// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dlarnv.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dorglq.dart';
import 'package:dart_lapack/src/dormlq.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';

import 'common.dart';

void dlqt03(
  final int M,
  final int N,
  final int K,
  final Matrix<double> AF_,
  final Matrix<double> C_,
  final Matrix<double> CC_,
  final Matrix<double> Q_,
  final int LDA,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<double> RESULT_,
) {
  final AF = AF_.having(ld: LDA);
  final C = C_.having(ld: LDA);
  final CC = CC_.having(ld: LDA);
  final Q = Q_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having(length: LWORK);
  final RWORK = RWORK_.having();
  final RESULT = RESULT_.having();
  const ONE = 1.0;
  const ROGUE = -1.0e+10;
  final ISEED = Array.fromList([1988, 1989, 1990, 1991]);
  final INFO = Box(0);

  final EPS = dlamch('Epsilon');

  // Copy the first k rows of the factorization to the array Q

  dlaset('Full', N, N, ROGUE, ROGUE, Q, LDA);
  dlacpy('Upper', K, N - 1, AF(1, 2), LDA, Q(1, 2), LDA);

  // Generate the n-by-n matrix Q

  srnamc.SRNAMT = 'DORGLQ';
  dorglq(N, N, K, Q, LDA, TAU, WORK, LWORK, INFO);

  for (var ISIDE = 1; ISIDE <= 2; ISIDE++) {
    final (SIDE, MC, NC) = ISIDE == 1 ? ('L', N, M) : ('R', M, N);

    // Generate MC by NC matrix C

    for (var J = 1; J <= NC; J++) {
      dlarnv(2, ISEED, MC, C(1, J).asArray());
    }
    final CNORM = switch (dlange('1', MC, NC, C, LDA, RWORK)) {
      0.0 => ONE,
      final i => i,
    };

    for (var ITRANS = 1; ITRANS <= 2; ITRANS++) {
      final TRANS = ITRANS == 1 ? 'N' : 'T';

      // Copy C

      dlacpy('Full', MC, NC, C, LDA, CC, LDA);

      // Apply Q or Q' to C

      srnamc.SRNAMT = 'DORMLQ';
      dormlq(SIDE, TRANS, MC, NC, K, AF, LDA, TAU, CC, LDA, WORK, LWORK, INFO);

      // Form explicit product and subtract

      if (lsame(SIDE, 'L')) {
        dgemm(TRANS, 'No transpose', MC, NC, MC, -ONE, Q, LDA, C, LDA, ONE, CC,
            LDA);
      } else {
        dgemm('No transpose', TRANS, MC, NC, NC, -ONE, C, LDA, Q, LDA, ONE, CC,
            LDA);
      }

      // Compute error in the difference

      final RESID = dlange('1', MC, NC, CC, LDA, RWORK);
      RESULT[(ISIDE - 1) * 2 + ITRANS] = RESID / (max(1, N) * CNORM * EPS);
    }
  }
}
