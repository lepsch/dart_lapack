// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';
import 'package:test/test.dart';

import '../matgen/dlatms.dart';
import '../test_driver.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'derrlq.dart';
import 'dget02.dart';
import 'dlarhs.dart';
import 'dlatb4.dart';
import 'dlqt01.dart';
import 'dlqt02.dart';
import 'dlqt03.dart';
import 'xlaenv.dart';

void dchklq(
  final Array<bool> DOTYPE_,
  final int NM,
  final Array<int> MVAL_,
  final int NN,
  final Array<int> NVAL_,
  final int NNB,
  final Array<int> NBVAL_,
  final Array<int> NXVAL_,
  final int NRHS,
  final double THRESH,
  final bool TSTERR,
  final int NMAX,
  final Array<double> A_,
  final Array<double> AF_,
  final Array<double> AQ_,
  final Array<double> AL_,
  final Array<double> AC_,
  final Array<double> B_,
  final Array<double> X_,
  final Array<double> XACT_,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final Array<double> RWORK_,
  final Nout NOUT,
  final TestDriver test,
) {
  final DOTYPE = DOTYPE_.having();
  final MVAL = MVAL_.having();
  final NVAL = NVAL_.having();
  final NBVAL = NBVAL_.having();
  final NXVAL = NXVAL_.having();
  final A = A_.having();
  final AF = AF_.having();
  final AQ = AQ_.having();
  final AL = AL_.having();
  final AC = AC_.having();
  final B = B_.having();
  final X = X_.having();
  final XACT = XACT_.having();
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const NTESTS = 7, NTYPES = 8;
  const ZERO = 0.0;
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}LQ';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  final ISEED = Array.fromList(ISEEDY);

  test.group('error exits', () {
    // Test the error exits
    if (TSTERR) derrlq(PATH, NOUT, test);
    test.tearDown(() {
      infoc.INFOT = 0;
    });
  });

  test.setUp(() {
    xlaenv(2, 2);
  });

  final LDA = NMAX;
  final LWORK = NMAX * max(NMAX, NRHS).toInt();

  // Do for each value of M in MVAL.
  for (final IM in 1.through(NM)) {
    final M = MVAL[IM];

    // Do for each value of N in NVAL.
    for (final IN in 1.through(NN)) {
      final N = NVAL[IN];
      final MINMN = min(M, N);
      for (final IMAT in 1.through(NTYPES)) {
        // Do the tests only if DOTYPE( IMAT ) is true.
        final skip = !DOTYPE[IMAT];

        test('DCHKLQ (IM=$IM IN=$IN IMAT=$IMAT)', () {
          final INFO = Box(0);

          // Set up parameters with DLATB4 and generate a test matrix
          // with DLATMS.
          final (:TYPE, :KL, :KU, :ANORM, :MODE, :COND, :DIST) =
              dlatb4(PATH, IMAT, M, N);

          srnamc.SRNAMT = 'DLATMS';
          dlatms(M, N, DIST, ISEED, TYPE, RWORK, MODE, COND, ANORM, KL, KU,
              'No packing', A.asMatrix(), LDA, WORK, INFO);

          // Check error code from DLATMS.
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            alaerh(PATH, 'DLATMS', INFO.value, 0, ' ', M, N, -1, -1, -1, IMAT,
                NFAIL, NERRS, NOUT);
            return;
          }

          // Set some values for K: the first value must be MINMN,
          // corresponding to the call of DLQT01; other values are
          // used in the calls of DLQT02, and must not exceed MINMN.
          final KVAL = Array.fromList([MINMN, 0, 1, MINMN ~/ 2]);
          final NK = switch (MINMN) {
            0 => 1,
            1 => 2,
            <= 3 => 3,
            _ => 4,
          };

          // Do for each value of K in KVAL
          for (var IK = 1; IK <= NK; IK++) {
            final K = KVAL[IK];

            // Do for each pair of values (NB,NX) in NBVAL and NXVAL.
            for (var INB = 1; INB <= NNB; INB++) {
              final NB = NBVAL[INB];
              xlaenv(1, NB);
              final NX = NXVAL[INB];
              xlaenv(3, NX);
              for (var I = 1; I <= NTESTS; I++) {
                RESULT[I] = ZERO;
              }
              var NT = 2;
              if (IK == 1) {
                // Test DGELQF
                dlqt01(M, N, A.asMatrix(), AF.asMatrix(), AQ.asMatrix(),
                    AL.asMatrix(), LDA, TAU, WORK, LWORK, RWORK, RESULT(1));
              } else if (M <= N) {
                // Test DORGLQ, using factorization returned by DLQT01
                dlqt02(M, N, K, A.asMatrix(), AF.asMatrix(), AQ.asMatrix(),
                    AL.asMatrix(), LDA, TAU, WORK, LWORK, RWORK, RESULT(1));
              } else {
                RESULT[1] = ZERO;
                RESULT[2] = ZERO;
              }
              if (M >= K) {
                // Test DORMLQ, using factorization returned by DLQT01
                dlqt03(M, N, K, AF.asMatrix(), AC.asMatrix(), AL.asMatrix(),
                    AQ.asMatrix(), LDA, TAU, WORK, LWORK, RWORK, RESULT(3));
                NT += 4;

                // If M<=N and K=M, call DGELS to solve a system
                // with NRHS right hand sides and compute the residual.
                if (K == M && INB == 1) {
                  // Generate a solution and set the right hand side.
                  srnamc.SRNAMT = 'DLARHS';
                  dlarhs(
                      PATH,
                      'New',
                      'Full',
                      'No transpose',
                      M,
                      N,
                      0,
                      0,
                      NRHS,
                      A.asMatrix(),
                      LDA,
                      XACT.asMatrix(),
                      LDA,
                      B.asMatrix(),
                      LDA,
                      ISEED,
                      INFO);

                  dlacpy('Full', M, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

                  // Reset AF to the original matrix. DGELS
                  // factors the matrix before solving the system.
                  dlacpy('Full', M, N, A.asMatrix(), LDA, AF.asMatrix(), LDA);

                  srnamc.SRNAMT = 'DGELS';
                  dgels('No transpose', M, N, NRHS, AF.asMatrix(), LDA,
                      X.asMatrix(), LDA, WORK, LWORK, INFO);

                  // Check error code from DGELS.
                  test.expect(INFO.value, 0);
                  if (INFO.value != 0) {
                    alaerh(PATH, 'DGELS', INFO.value, 0, 'N', M, N, NRHS, -1,
                        NB, IMAT, NFAIL, NERRS, NOUT);
                  }

                  dget02('No transpose', M, N, NRHS, A.asMatrix(), LDA,
                      X.asMatrix(), LDA, B.asMatrix(), LDA, RWORK, RESULT(7));
                  NT++;
                } else {
                  RESULT[7] = ZERO;
                }
              } else {
                RESULT[3] = ZERO;
                RESULT[4] = ZERO;
                RESULT[5] = ZERO;
                RESULT[6] = ZERO;
              }

              // Print information about the tests that did not
              // pass the threshold.

              for (var I = 1; I <= NT; I++) {
                final reason =
                    ' M=${M.i5}, N=${N.i5}, K=${K.i5}, NB=${NB.i4}, NX=${NX.i5}, type ${IMAT.i2}, test(${I.i2})=${RESULT[I].g12_5}';
                test.expect(RESULT[I], lessThan(THRESH), reason: reason);
                if (RESULT[I] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                  NOUT.println(reason);
                  NFAIL++;
                }
              }
              NRUN += NT;
            }
          }
        }, skip: skip);
      }
    }
  }

  // Print a summary of the results.
  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
