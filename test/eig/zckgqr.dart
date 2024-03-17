import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../matgen/zlatms.dart';
import 'alahdg.dart';
import 'alareq.dart';
import 'alasum.dart';
import 'dlatb9.dart';
import 'zgqrts.dart';
import 'zgrqts.dart';

Future<void> zckgqr(
  final int NM,
  final Array<int> MVAL_,
  final int NP,
  final Array<int> PVAL_,
  final int NN,
  final Array<int> NVAL_,
  final int NMATS,
  final Array<int> ISEED_,
  final double THRESH,
  final int NMAX,
  final Array<Complex> A_,
  final Array<Complex> AF_,
  final Array<Complex> AQ_,
  final Array<Complex> AR_,
  final Array<Complex> TAUA_,
  final Array<Complex> B_,
  final Array<Complex> BF_,
  final Array<Complex> BZ_,
  final Array<Complex> BT_,
  final Array<Complex> BWK_,
  final Array<Complex> TAUB_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Nin NIN,
  final Nout NOUT,
  final Box<int> INFO,
) async {
  final A = A_.having();
  final AF = AF_.having();
  final AQ = AQ_.having();
  final AR = AR_.having();
  final TAUA = TAUA_.having();
  final B = B_.having();
  final BF = BF_.having();
  final BZ = BZ_.having();
  final BT = BT_.having();
  final BWK = BWK_.having();
  final TAUB = TAUB_.having();
  final MVAL = MVAL_.having();
  final PVAL = PVAL_.having();
  final NVAL = NVAL_.having();
  final ISEED = ISEED_.having(length: 4);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NTESTS = 7;
  const NTYPES = 8;
  bool FIRSTT;
  String PATH;
  int I, IM, IMAT, IN, IP, LDA, LDB, LWORK, M, N, NT, P;
  final DOTYPE = Array<bool>(NTYPES);
  final RESULT = Array<double>(NTESTS);
  final DISTA = Box(''), DISTB = Box(''), TYPE = Box('');
  final ANORM = Box(0.0),
      BNORM = Box(0.0),
      CNDNMA = Box(0.0),
      CNDNMB = Box(0.0);
  final IINFO = Box(0),
      KLA = Box(0),
      KLB = Box(0),
      KUA = Box(0),
      KUB = Box(0),
      MODEA = Box(0),
      MODEB = Box(0);

  // Initialize constants.

  PATH = 'GQR';
  INFO.value = 0;
  var NRUN = 0;
  var NFAIL = 0;
  FIRSTT = true;
  await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);
  LDA = NMAX;
  LDB = NMAX;
  LWORK = NMAX * NMAX;

  // Do for each value of M in MVAL.

  for (IM = 1; IM <= NM; IM++) {
    M = MVAL[IM];

    // Do for each value of P in PVAL.

    for (IP = 1; IP <= NP; IP++) {
      P = PVAL[IP];

      // Do for each value of N in NVAL.

      for (IN = 1; IN <= NN; IN++) {
        N = NVAL[IN];

        for (IMAT = 1; IMAT <= NTYPES; IMAT++) {
          // Do the tests only if DOTYPE( IMAT ) is true.

          if (!DOTYPE[IMAT]) continue;

          // Test ZGGRQF

          // Set up parameters with DLATB9 and generate test
          // matrices A and B with ZLATMS.

          dlatb9('GRQ', IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, ANORM, BNORM,
              MODEA, MODEB, CNDNMA, CNDNMB, DISTA, DISTB);

          zlatms(
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
              A.asMatrix(),
              LDA,
              WORK,
              IINFO);
          if (IINFO.value != 0) {
            _print9999(NOUT, IINFO.value);
            INFO.value = (IINFO.value).abs();
            continue;
          }

          zlatms(
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
              B.asMatrix(),
              LDB,
              WORK,
              IINFO);
          if (IINFO.value != 0) {
            _print9999(NOUT, IINFO.value);
            INFO.value = (IINFO.value).abs();
            continue;
          }

          NT = 4;

          zgrqts(
              M,
              P,
              N,
              A.asMatrix(),
              AF.asMatrix(),
              AQ.asMatrix(),
              AR.asMatrix(),
              LDA,
              TAUA,
              B.asMatrix(),
              BF.asMatrix(),
              BZ.asMatrix(),
              BT.asMatrix(),
              BWK.asMatrix(),
              LDB,
              TAUB,
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
                alahdg(NOUT, 'GRQ');
              }
              NOUT.println(
                  ' M=${M.i4} P=${P.i4}, N=${N.i4}, type ${IMAT.i2}, test ${I.i2}, ratio=${RESULT[I].g13_6}');
              NFAIL++;
            }
          }
          NRUN += NT;

          // Test ZGGQRF

          // Set up parameters with DLATB9 and generate test
          // matrices A and B with ZLATMS.

          dlatb9('GQR', IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, ANORM, BNORM,
              MODEA, MODEB, CNDNMA, CNDNMB, DISTA, DISTB);

          zlatms(
              N,
              M,
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
              A.asMatrix(),
              LDA,
              WORK,
              IINFO);
          if (IINFO.value != 0) {
            _print9999(NOUT, IINFO.value);
            INFO.value = (IINFO.value).abs();
            continue;
          }

          zlatms(
              N,
              P,
              DISTB.value,
              ISEED,
              TYPE.value,
              RWORK,
              MODEA.value,
              CNDNMA.value,
              BNORM.value,
              KLB.value,
              KUB.value,
              'No packing',
              B.asMatrix(),
              LDB,
              WORK,
              IINFO);
          if (IINFO.value != 0) {
            _print9999(NOUT, IINFO.value);
            INFO.value = (IINFO.value).abs();
            continue;
          }

          NT = 4;

          zgqrts(
              N,
              M,
              P,
              A.asMatrix(),
              AF.asMatrix(),
              AQ.asMatrix(),
              AR.asMatrix(),
              LDA,
              TAUA,
              B.asMatrix(),
              BF.asMatrix(),
              BZ.asMatrix(),
              BT.asMatrix(),
              BWK.asMatrix(),
              LDB,
              TAUB,
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
                  ' N=${N.i4} M=${M.i4}, P=${P.i4}, type ${IMAT.i2}, test ${I.i2}, ratio=${RESULT[I].g13_6}');
              NFAIL++;
            }
          }
          NRUN += NT;
        }
      }
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, 0);
}

void _print9999(Nout nout, int info) {
  nout.println(' ZLATMS in ZCKGQR:    INFO = ${info.i5}');
}
