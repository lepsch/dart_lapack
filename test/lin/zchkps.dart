// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zpstrf.dart';

import '../matgen/zlatmt.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'xlaenv.dart';
import 'zerrps.dart';
import 'zlatb5.dart';
import 'zpst01.dart';

void zchkps(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NNB,
  final Array<int> NBVAL_,
  final int NRANK,
  final Array<int> RANKVAL_,
  final double THRESH,
  final bool TSTERR,
  final int NMAX,
  final Array<Complex> A_,
  final Array<Complex> AFAC_,
  final Array<Complex> PERM_,
  final Array<int> PIV_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Nout NOUT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DOTYPE = DOTYPE_.having();
  final NVAL = NVAL_.having();
  final NBVAL = NBVAL_.having();
  final RANKVAL = RANKVAL_.having();
  final A = A_.having();
  final AFAC = AFAC_.having();
  final PERM = PERM_.having();
  final PIV = PIV_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();

  const ONE = 1.0;
  const NTYPES = 9;
  final ISEED = Array<int>(4);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Zomplex Precision'[0]}PS';
  var NRUN = 0;
  var NFAIL = 0;
  var NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) zerrps(PATH, NOUT);
  infoc.INFOT = 0;

  // Do for each value of N in NVAL

  for (var IN = 1; IN <= NN; IN++) {
    final N = NVAL[IN];
    final LDA = max(N, 1);
    final NIMAT = N <= 0 ? 1 : NTYPES;

    for (var IMAT = 1; IMAT <= NIMAT; IMAT++) {
      // Do the tests only if DOTYPE( IMAT ) is true.

      if (!DOTYPE[IMAT]) continue;

      // Do for each value of RANK in RANKVAL

      for (var IRANK = 1; IRANK <= NRANK; IRANK++) {
        // Only repeat test 3 to 5 for different ranks
        // Other tests use full rank

        if ((IMAT < 3 || IMAT > 5) && IRANK > 1) continue;

        final RANK = ((N * RANKVAL[IRANK]) / 100.0).ceil();

        // Do first for UPLO = 'U', then for UPLO = 'L'

        for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
          final UPLO = UPLOS[IUPLO - 1];

          // Set up parameters with ZLATB5 and generate a test matrix
          // with ZLATMT.

          final (:TYPE, :KL, :KU, :ANORM, :MODE, COND: CNDNUM, :DIST) =
              zlatb5(PATH, IMAT, N);

          srnamc.SRNAMT = 'ZLATMT';
          zlatmt(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, RANK, KL,
              KU, UPLO, A.asMatrix(), LDA, WORK, INFO);

          // Check error code from ZLATMT.

          if (INFO.value != 0) {
            alaerh(PATH, 'ZLATMT', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
                NFAIL, NERRS, NOUT);
            continue;
          }

          // Do for each value of NB in NBVAL

          for (var INB = 1; INB <= NNB; INB++) {
            final NB = NBVAL[INB];
            xlaenv(1, NB);

            // Compute the pivoted L*L' or U'*U factorization
            // of the matrix.

            zlacpy(UPLO, N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);
            srnamc.SRNAMT = 'ZPSTRF';

            // Use default tolerance

            final TOL = -ONE;
            final COMPRANK = Box(0);
            zpstrf(
                UPLO, N, AFAC.asMatrix(), LDA, PIV, COMPRANK, TOL, RWORK, INFO);

            // Check error code from ZPSTRF.

            const IZERO = 0;
            if ((INFO.value < IZERO) ||
                (INFO.value != IZERO && RANK == N) ||
                (INFO.value <= IZERO && RANK < N)) {
              alaerh(PATH, 'ZPSTRF', INFO.value, IZERO, UPLO, N, N, -1, -1, NB,
                  IMAT, NFAIL, NERRS, NOUT);
              continue;
            }

            // Skip the test if INFO is not 0.

            if (INFO.value != 0) continue;

            // Reconstruct matrix from factors and compute residual.

            // PERM holds permuted L*L^T or U^T*U

            final RESULT = Box(0.0);
            zpst01(UPLO, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA,
                PERM.asMatrix(), LDA, PIV, RWORK, RESULT, COMPRANK.value);

            // Print information about the tests that did not pass
            // the threshold or where computed rank was not RANK.

            if (N == 0) COMPRANK.value = 0;
            final RANKDIFF = RANK - COMPRANK.value;
            if (RESULT.value >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.println(
                  ' UPLO = \'${UPLO.a1}\', N =${N.i5}, RANK =${RANK.i3}, Diff =${RANKDIFF.i5}, NB =${NB.i4}, type ${IMAT.i2}, Ratio =${RESULT.value.g12_5}');
              NFAIL++;
            }
            NRUN++;
          }
        }
      }
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
