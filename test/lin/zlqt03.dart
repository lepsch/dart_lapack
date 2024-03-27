import 'dart:math';

import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlarnv.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zunglq.dart';
import 'package:lapack/src/zunmlq.dart';

import 'common.dart';

void zlqt03(
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
  final AF = AF_.having(ld: LDA);
  final C = C_.having(ld: LDA);
  final CC = CC_.having(ld: LDA);
  final Q = Q_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having(length: LWORK);
  final RWORK = RWORK_.having();
  final RESULT = RESULT_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0;
  const ROGUE = Complex(-1.0e+10, -1.0e+10);
  final ISEED = Array.fromList([1988, 1989, 1990, 1991]);
  final INFO = Box(0);

  final EPS = dlamch('Epsilon');

  // Copy the first k rows of the factorization to the array Q

  zlaset('Full', N, N, ROGUE, ROGUE, Q, LDA);
  zlacpy('Upper', K, N - 1, AF(1, 2), LDA, Q(1, 2), LDA);

  // Generate the n-by-n matrix Q

  srnamc.SRNAMT = 'ZUNGLQ';
  zunglq(N, N, K, Q, LDA, TAU, WORK, LWORK, INFO);

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

      srnamc.SRNAMT = 'ZUNMLQ';
      zunmlq(SIDE, TRANS, MC, NC, K, AF, LDA, TAU, CC, LDA, WORK, LWORK, INFO);

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
