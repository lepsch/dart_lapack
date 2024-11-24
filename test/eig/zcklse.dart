// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';

import '../lin/zlarhs.dart';
import '../matgen/zlatms.dart';
import 'alahdg.dart';
import 'alareq.dart';
import 'alasum.dart';
import 'dlatb9.dart';
import 'zlsets.dart';

Future<void> zcklse(
  final int NN,
  final Array<int> MVAL_,
  final Array<int> PVAL_,
  final Array<int> NVAL_,
  final int NMATS,
  final Array<int> ISEED_,
  final double THRESH,
  final int NMAX,
  final Array<Complex> A_,
  final Array<Complex> AF_,
  final Array<Complex> B_,
  final Array<Complex> BF_,
  final Array<Complex> X_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Nin NIN,
  final Nout NOUT,
  final Box<int> INFO,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ISEED = ISEED_.having(length: 4);
  final A = A_.having();
  final AF = AF_.having();
  final B = B_.having();
  final BF = BF_.having();
  final X = X_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final MVAL = MVAL_.having();
  final PVAL = PVAL_.having();
  final NVAL = NVAL_.having();
  const NTESTS = 7;
  const NTYPES = 8;
  bool FIRSTT;
  String PATH;
  int I, IK, IMAT, LDA, LDB, LWORK, M, N, NT, P;
  final DOTYPE = Array<bool>(NTYPES);
  final RESULT = Array<double>(NTESTS);
  final IINFO = Box(0);

  // Initialize constants and the random number seed.

  PATH = 'LSE';
  INFO.value = 0;
  var NRUN = 0;
  var NFAIL = 0;
  FIRSTT = true;
  await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);
  LDA = NMAX;
  LDB = NMAX;
  LWORK = NMAX * NMAX;

  // Check for valid input values.

  for (IK = 1; IK <= NN; IK++) {
    M = MVAL[IK];
    P = PVAL[IK];
    N = NVAL[IK];
    if (P > N || N > M + P) {
      if (FIRSTT) {
        NOUT.println();
        FIRSTT = false;
      }
      NOUT.println(
          ' *** Invalid input  for LSE:  M = ${M.i6}, P = ${P.i6}, N = ${N.i6};\n     must satisfy P <= N <= P+M  (this set of values will be skipped)');
    }
  }
  FIRSTT = true;

  // Do for each value of M in MVAL.

  for (IK = 1; IK <= NN; IK++) {
    M = MVAL[IK];
    P = PVAL[IK];
    N = NVAL[IK];
    if (P > N || N > M + P) continue;

    for (IMAT = 1; IMAT <= NTYPES; IMAT++) {
      // Do the tests only if DOTYPE[ IMAT ] is true.

      if (!DOTYPE[IMAT]) continue;

      // Set up parameters with DLATB9 and generate test
      // matrices A and B with ZLATMS.

      final (
        :TYPE,
        :KLA,
        :KUA,
        :KLB,
        :KUB,
        :ANORM,
        :BNORM,
        :MODEA,
        :MODEB,
        :CNDNMA,
        :CNDNMB,
        :DISTA,
        :DISTB
      ) = dlatb9(PATH, IMAT, M, P, N);

      zlatms(M, N, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, ANORM, KLA, KUA,
          'No packing', A.asMatrix(), LDA, WORK, IINFO);
      if (IINFO.value != 0) {
        _print9999(NOUT, IINFO.value);
        INFO.value = (IINFO.value).abs();
        continue;
      }

      zlatms(P, N, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB, BNORM, KLB, KUB,
          'No packing', B.asMatrix(), LDB, WORK, IINFO);
      if (IINFO.value != 0) {
        _print9999(NOUT, IINFO.value);
        INFO.value = (IINFO.value).abs();
        continue;
      }

      // Generate the right-hand sides C and D for the LSE.

      zlarhs(
          'ZGE',
          'New solution',
          'Upper',
          'N',
          M,
          N,
          max(M - 1, 0),
          max(N - 1, 0),
          1,
          A.asMatrix(),
          LDA,
          X(4 * NMAX + 1).asMatrix(),
          max(N, 1),
          X.asMatrix(),
          max(M, 1),
          ISEED,
          IINFO);

      zlarhs(
          'ZGE',
          'Computed',
          'Upper',
          'N',
          P,
          N,
          max(P - 1, 0),
          max(N - 1, 0),
          1,
          B.asMatrix(),
          LDB,
          X(4 * NMAX + 1).asMatrix(),
          max(N, 1),
          X(2 * NMAX + 1).asMatrix(),
          max(P, 1),
          ISEED,
          IINFO);

      NT = 2;

      zlsets(
          M,
          P,
          N,
          A.asMatrix(),
          AF.asMatrix(),
          LDA,
          B.asMatrix(),
          BF.asMatrix(),
          LDB,
          X,
          X(NMAX + 1),
          X(2 * NMAX + 1),
          X(3 * NMAX + 1),
          X(4 * NMAX + 1),
          WORK,
          LWORK,
          RWORK,
          RESULT(1));

      // Print information about the tests that did not
      // pass the threshold.

      for (I = 1; I <= NT; I++) {
        if (RESULT[I] >= THRESH) {
          if (NFAIL == 0 && FIRSTT) {
            FIRSTT = false;
            alahdg(NOUT, PATH);
          }
          NOUT.println(
              ' M=${M.i4} P=${P.i4}, N=${N.i4}, type ${IMAT.i2}, test ${I.i2}, ratio=${RESULT[I].g13_6}');
          NFAIL++;
        }
      }
      NRUN += NT;
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, 0);
}

void _print9999(Nout nout, int info) {
  nout.println(' ZLATMS in ZCKLSE   INFO = ${info.i5}');
}
