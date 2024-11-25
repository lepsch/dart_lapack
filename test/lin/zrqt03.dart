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
import 'package:dart_lapack/src/zungrq.dart';
import 'package:dart_lapack/src/zunmrq.dart';

import 'common.dart';

void zrqt03(
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
  final AF = AF_.having();
  final C = C_.having();
  final CC = CC_.having();
  final Q = Q_.having();
  final TAU = TAU_.having();
  final WORK = WORK_.having(length: LWORK);
  final RWORK = RWORK_.having();
  final RESULT = RESULT_.having();
  const ZERO = 0.0, ONE = 1.0;
  const ROGUE = Complex(-1.0e+10, -1.0e+10);
  final INFO = Box(0);
  final ISEED = Array.fromList([1988, 1989, 1990, 1991]);

  final EPS = dlamch('Epsilon');
  final MINMN = min(M, N);

  // Quick return if possible

  if (MINMN == 0) {
    RESULT[1] = ZERO;
    RESULT[2] = ZERO;
    RESULT[3] = ZERO;
    RESULT[4] = ZERO;
    return;
  }

  // Copy the last k rows of the factorization to the array Q

  zlaset('Full', N, N, ROGUE, ROGUE, Q, LDA);
  if (K > 0 && N > K) {
    zlacpy('Full', K, N - K, AF(M - K + 1, 1), LDA, Q(N - K + 1, 1), LDA);
  }
  if (K > 1) {
    zlacpy('Lower', K - 1, K - 1, AF(M - K + 2, N - K + 1), LDA,
        Q(N - K + 2, N - K + 1), LDA);
  }

  // Generate the n-by-n matrix Q

  srnamc.SRNAMT = 'ZUNGRQ';
  zungrq(N, N, K, Q, LDA, TAU(MINMN - K + 1), WORK, LWORK, INFO);

  for (var ISIDE = 1; ISIDE <= 2; ISIDE++) {
    final (SIDE, MC, NC) = ISIDE == 1 ? ('L', N, M) : ('R', M, N);

    // Generate MC by NC matrix C

    for (var J = 1; J <= NC; J++) {
      zlarnv(2, ISEED, MC, C(1, J).asArray());
    }
    final CNORM = switch (zlange('1', MC, NC, C, LDA, RWORK)) {
      0 => ONE,
      final d => d,
    };

    for (var ITRANS = 1; ITRANS <= 2; ITRANS++) {
      final TRANS = ITRANS == 1 ? 'N' : 'C';

      // Copy C

      zlacpy('Full', MC, NC, C, LDA, CC, LDA);

      // Apply Q or Q' to C

      srnamc.SRNAMT = 'ZUNMRQ';
      if (K > 0) {
        zunmrq(SIDE, TRANS, MC, NC, K, AF(M - K + 1, 1), LDA,
            TAU(MINMN - K + 1), CC, LDA, WORK, LWORK, INFO);
      }

      // Form explicit product and subtract

      if (lsame(SIDE, 'L')) {
        zgemm(TRANS, 'No transpose', MC, NC, MC, Complex(-ONE), Q, LDA, C, LDA,
            Complex.one, CC, LDA);
      } else {
        zgemm('No transpose', TRANS, MC, NC, NC, Complex(-ONE), C, LDA, Q, LDA,
            Complex.one, CC, LDA);
      }

      // Compute error in the difference

      final RESID = zlange('1', MC, NC, CC, LDA, RWORK);
      RESULT[(ISIDE - 1) * 2 + ITRANS] = RESID / (max(1, N) * CNORM * EPS);
    }
  }
}
