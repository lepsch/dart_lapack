import 'package:lapack/src/box.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../lin/alasum.dart';
import '../matgen/dlatms.dart';
import 'alahdg.dart';
import 'alareq.dart';
import 'dgsvts3.dart';
import 'dlatb9.dart';

Future<void> dckgsv(
  final int NM,
  final Array<int> MVAL,
  final Array<int> PVAL,
  final Array<int> NVAL,
  final int NMATS,
  final Array<int> ISEED,
  final double THRESH,
  final int NMAX,
  final Array<double> A,
  final Array<double> AF,
  final Array<double> B,
  final Array<double> BF,
  final Array<double> U,
  final Array<double> V,
  final Array<double> Q,
  final Array<double> ALPHA,
  final Array<double> BETA,
  final Array<double> R,
  final Array<int> IWORK,
  final Array<double> WORK,
  final Array<double> RWORK,
  final Nin NIN,
  final Nout NOUT,
  final Box<int> INFO,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NTESTS = 12;
  const NTYPES = 8;
  bool FIRSTT;
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
  const PATH = 'GSV';
  final DISTA = Box(''), DISTB = Box(''), TYPE = Box('');
  final KLA = Box(0),
      KLB = Box(0),
      KUA = Box(0),
      KUB = Box(0),
      MODEA = Box(0),
      MODEB = Box(0),
      IINFO = Box(0);
  final ANORM = Box(0.0),
      BNORM = Box(0.0),
      CNDNMA = Box(0.0),
      CNDNMB = Box(0.0);

  // Initialize constants and the random number seed.

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
      // matrices A and B with DLATMS.

      dlatb9(
        PATH,
        IMAT,
        M,
        P,
        N,
        TYPE,
        KLA,
        KUA,
        KLB,
        KUB,
        ANORM,
        BNORM,
        MODEA,
        MODEB,
        CNDNMA,
        CNDNMB,
        DISTA,
        DISTB,
      );

      // Generate M by N matrix A

      dlatms(
        M,
        N,
        DISTA.value,
        ISEED,
        TYPE.value,
        RWORK,
        MODEA.value,
        CNDNMA.value,
        ANORM.value,
        KLA.value,
        KUA.value,
        'No packing',
        A.asMatrix(LDA),
        LDA,
        WORK,
        IINFO,
      );
      if (IINFO.value != 0) {
        print9999(NOUT, IINFO.value);
        INFO.value = (IINFO.value).abs();
        continue;
      }

      dlatms(
        P,
        N,
        DISTB.value,
        ISEED,
        TYPE.value,
        RWORK,
        MODEB.value,
        CNDNMB.value,
        BNORM.value,
        KLB.value,
        KUB.value,
        'No packing',
        B.asMatrix(LDB),
        LDB,
        WORK,
        IINFO,
      );
      if (IINFO.value != 0) {
        print9999(NOUT, IINFO.value);
        INFO.value = (IINFO.value).abs();
        continue;
      }

      NT = 6;

      dgsvts3(
        M,
        P,
        N,
        A,
        AF,
        LDA,
        B,
        BF,
        LDB,
        U,
        LDU,
        V,
        LDV,
        Q,
        LDQ,
        ALPHA,
        BETA,
        R,
        LDR,
        IWORK,
        WORK,
        LWORK,
        RWORK,
        RESULT,
      );

      // Print information about the tests that did not
      // pass the threshold.

      for (I = 1; I <= NT; I++) {
        if (RESULT[I] >= THRESH) {
          if (NFAIL == 0 && FIRSTT) {
            FIRSTT = false;
            alahdg(NOUT, PATH);
          }
          NOUT.println(
            ' M=${M.i4} P=${P.i4}, N=${N.i4}, type ${IMAT.i2}, test ${I.i2}, ratio=${RESULT[I].g13_6}',
          );
          NFAIL = NFAIL + 1;
        }
      }
      NRUN = NRUN + NT;
    }
  }

  // Print a summary of the results.
  alasum(PATH, NOUT, NFAIL, NRUN, 0);
}

void print9999(final Nout NOUT, final int info) {
  NOUT.println(' DLATMS in DCKGSV   INFO.value = ${info.i5}');
}
