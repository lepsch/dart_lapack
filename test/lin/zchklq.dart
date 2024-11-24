// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/zgels.dart';
import 'package:dart_lapack/src/zlacpy.dart';

import '../matgen/zlatms.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'xlaenv.dart';
import 'zerrlq.dart';
import 'zget02.dart';
import 'zlarhs.dart';
import 'zlatb4.dart';
import 'zlqt01.dart';
import 'zlqt02.dart';
import 'zlqt03.dart';

void zchklq(
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
  final Array<Complex> A_,
  final Array<Complex> AF_,
  final Array<Complex> AQ_,
  final Array<Complex> AL_,
  final Array<Complex> AC_,
  final Array<Complex> B_,
  final Array<Complex> X_,
  final Array<Complex> XACT_,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Nout NOUT,
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

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NTESTS = 7, NTYPES = 8;
  const ZERO = 0.0;
  final ISEED = Array<int>(4), KVAL = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Zomplex precision'[0]}LQ';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) zerrlq(PATH, NOUT);
  infoc.INFOT = 0;
  xlaenv(2, 2);

  final LDA = NMAX;
  final LWORK = NMAX * max(NMAX, NRHS).toInt();

  // Do for each value of M in MVAL.

  for (var IM = 1; IM <= NM; IM++) {
    final M = MVAL[IM];

    // Do for each value of N in NVAL.

    for (var IN = 1; IN <= NN; IN++) {
      final N = NVAL[IN];
      final MINMN = min(M, N);
      for (var IMAT = 1; IMAT <= NTYPES; IMAT++) {
        // Do the tests only if DOTYPE( IMAT ) is true.

        if (!DOTYPE[IMAT]) continue;

        // Set up parameters with ZLATB4 and generate a test matrix
        // with ZLATMS.

        final (:TYPE, :KL, :KU, :ANORM, :MODE, :CNDNUM, :DIST) =
            zlatb4(PATH, IMAT, M, N);

        srnamc.SRNAMT = 'ZLATMS';
        zlatms(M, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
            'No packing', A.asMatrix(), LDA, WORK, INFO);

        // Check error code from ZLATMS.

        if (INFO.value != 0) {
          alaerh(PATH, 'ZLATMS', INFO.value, 0, ' ', M, N, -1, -1, -1, IMAT,
              NFAIL, NERRS, NOUT);
          continue;
        }

        // Set some values for K: the first value must be MINMN,
        // corresponding to the call of ZLQT01; other values are
        // used in the calls of ZLQT02, and must not exceed MINMN.

        KVAL[1] = MINMN;
        KVAL[2] = 0;
        KVAL[3] = 1;
        KVAL[4] = MINMN ~/ 2;
        final int NK;
        if (MINMN == 0) {
          NK = 1;
        } else if (MINMN == 1) {
          NK = 2;
        } else if (MINMN <= 3) {
          NK = 3;
        } else {
          NK = 4;
        }

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
              // Test ZGELQF

              zlqt01(M, N, A.asMatrix(), AF.asMatrix(), AQ.asMatrix(),
                  AL.asMatrix(), LDA, TAU, WORK, LWORK, RWORK, RESULT(1));
            } else if (M <= N) {
              // Test ZUNGLQ, using factorization
              // returned by ZLQT01

              zlqt02(M, N, K, A.asMatrix(), AF.asMatrix(), AQ.asMatrix(),
                  AL.asMatrix(), LDA, TAU, WORK, LWORK, RWORK, RESULT(1));
            }
            if (M >= K) {
              // Test ZUNMLQ, using factorization returned
              // by ZLQT01

              zlqt03(M, N, K, AF.asMatrix(), AC.asMatrix(), AL.asMatrix(),
                  AQ.asMatrix(), LDA, TAU, WORK, LWORK, RWORK, RESULT(3));
              NT += 4;

              // If M<=N and K=M, call ZGELS to solve a system
              // with NRHS right hand sides and compute the
              // residual.

              if (K == M && INB == 1) {
                // Generate a solution and set the right
                // hand side.

                srnamc.SRNAMT = 'ZLARHS';
                zlarhs(
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

                zlacpy('Full', M, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

                // Reset AF to the original matrix. ZGELS
                // factors the matrix before solving the system.

                zlacpy('Full', M, N, A.asMatrix(), LDA, AF.asMatrix(), LDA);

                srnamc.SRNAMT = 'ZGELS';
                zgels('No transpose', M, N, NRHS, AF.asMatrix(), LDA,
                    X.asMatrix(), LDA, WORK, LWORK, INFO);

                // Check error code from ZGELS.

                if (INFO.value != 0) {
                  alaerh(PATH, 'ZGELS', INFO.value, 0, 'N', M, N, NRHS, -1, NB,
                      IMAT, NFAIL, NERRS, NOUT);
                }

                zget02('No transpose', M, N, NRHS, A.asMatrix(), LDA,
                    X.asMatrix(), LDA, B.asMatrix(), LDA, RWORK, RESULT(7));
                NT++;
              }
            }

            // Print information about the tests that did not
            // pass the threshold.

            for (var I = 1; I <= NT; I++) {
              if (RESULT[I] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                NOUT.println(
                    ' M=${M.i5}, N=${N.i5}, K=${K.i5}, NB=${NB.i4}, NX=${NX.i5}, type ${IMAT.i2}, test(${I.i2})=${RESULT[I].g12_5}');
                NFAIL++;
              }
            }
            NRUN += NT;
          }
        }
      }
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
