// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';

import '../matgen/zlatms.dart';
import 'alahdg.dart';
import 'alareq.dart';
import 'alasum.dart';
import 'dlatb9.dart';
import 'zgsvts3.dart';

Future<void> zckgsv(
  final int NM,
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
  final Array<Complex> U_,
  final Array<Complex> V_,
  final Array<Complex> Q_,
  final Array<double> ALPHA_,
  final Array<double> BETA_,
  final Array<Complex> R_,
  final Array<int> IWORK_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Nin NIN,
  final Nout NOUT,
  final Box<int> INFO,
) async {
  final A = A_.having();
  final AF = AF_.having();
  final B = B_.having();
  final BF = BF_.having();
  final U = U_.having();
  final V = V_.having();
  final Q = Q_.having();
  final R = R_.having();
  final MVAL = MVAL_.having();
  final PVAL = PVAL_.having();
  final NVAL = NVAL_.having();
  final ALPHA = ALPHA_.having();
  final BETA = BETA_.having();
  final ISEED = ISEED_.having(length: 4);
  final IWORK = IWORK_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();

  const NTESTS = 12;
  const NTYPES = 8;
  bool FIRSTT;
  String PATH;
  int I,
      IM,
      IMAT,
      LDA,
      LDB,
      LDQ,
      LDR,
      LDU,
      LDV,
      LWORK,
      M,
      N,
      NFAIL,
      NRUN,
      NT,
      P;
  final DOTYPE = Array<bool>(NTYPES);
  final RESULT = Array<double>(NTESTS);
  final IINFO = Box(0);

  // Initialize constants and the random number seed.

  PATH = 'GSV';
  INFO.value = 0;
  NRUN = 0;
  NFAIL = 0;
  FIRSTT = true;
  await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);
  LDA = NMAX;
  LDB = NMAX;
  LDU = NMAX;
  LDV = NMAX;
  LDQ = NMAX;
  LDR = NMAX;
  LWORK = NMAX * NMAX;

  // Do for each value of M in MVAL.

  for (IM = 1; IM <= NM; IM++) {
    M = MVAL[IM];
    P = PVAL[IM];
    N = NVAL[IM];

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

      // Generate M by N matrix A

      zlatms(M, N, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, ANORM, KLA, KUA,
          'No packing', A.asMatrix(), LDA, WORK, IINFO);
      if (IINFO.value != 0) {
        _print9999(NOUT, IINFO.value);
        INFO.value = (IINFO.value).abs();
        continue;
      }

      // Generate P by N matrix B

      zlatms(P, N, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB, BNORM, KLB, KUB,
          'No packing', B.asMatrix(), LDB, WORK, IINFO);
      if (IINFO.value != 0) {
        _print9999(NOUT, IINFO.value);
        INFO.value = (IINFO.value).abs();
        continue;
      }

      NT = 6;

      zgsvts3(
          M,
          P,
          N,
          A.asMatrix(),
          AF.asMatrix(),
          LDA,
          B.asMatrix(),
          BF.asMatrix(),
          LDB,
          U.asMatrix(),
          LDU,
          V.asMatrix(),
          LDV,
          Q.asMatrix(),
          LDQ,
          ALPHA,
          BETA,
          R.asMatrix(),
          LDR,
          IWORK,
          WORK,
          LWORK,
          RWORK,
          RESULT);

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
  nout.println(' ZLATMS in ZCKGSV   INFO = ${info.i5}');
}
