import 'package:lapack/src/box.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../lin/alasum.dart';
import '../matgen/dlatms.dart';
import 'alahdg.dart';
import 'alareq.dart';
import 'dgqrts.dart';
import 'dgrqts.dart';
import 'dlatb9.dart';

Future<void> dckgqr(
  final int NM,
  final Array<int> MVAL,
  final int NP,
  final Array<int> PVAL,
  final int NN,
  final Array<int> NVAL,
  final int NMATS,
  final Array<int> ISEED,
  final double THRESH,
  final int NMAX,
  final Array<double> A,
  final Array<double> AF,
  final Array<double> AQ,
  final Array<double> AR,
  final Array<double> TAUA,
  final Array<double> B,
  final Array<double> BF,
  final Array<double> BZ,
  final Array<double> BT,
  final Array<double> BWK,
  final Array<double> TAUB,
  final Array<double> WORK,
  final Array<double> RWORK,
  final Nin NIN,
  final Nout NOUT,
  final Box<int> INFO,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NTESTS = 7;
  const NTYPES = 8;
  bool FIRSTT;
  final DISTA = Box(''), DISTB = Box(''), TYPE = Box('');
  int I, IM, IMAT, IN, IP, LDA, LDB, LWORK, M, N, NFAIL, NRUN, NT, P;
  final ANORM = Box(0.0),
      BNORM = Box(0.0),
      CNDNMA = Box(0.0),
      CNDNMB = Box(0.0);
  final DOTYPE = Array<bool>(NTYPES);
  final RESULT = Array<double>(NTESTS);
  final IINFO = Box(0),
      KLA = Box(0),
      KLB = Box(0),
      KUA = Box(0),
      KUB = Box(0),
      MODEA = Box(0),
      MODEB = Box(0);
  const PATH = 'GQR';

  // Initialize constants.

  INFO.value = 0;
  NRUN = 0;
  NFAIL = 0;
  FIRSTT = true;
  await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);
  LDA = NMAX;
  LDB = NMAX;
  LWORK = NMAX * NMAX;

  // Do for each value of M in MVAL.

  for (IM = 1; IM <= NM; IM++) {
    // 60
    M = MVAL[IM];

    // Do for each value of P in PVAL.

    for (IP = 1; IP <= NP; IP++) {
      // 50
      P = PVAL[IP];

      // Do for each value of N in NVAL.

      for (IN = 1; IN <= NN; IN++) {
        // 40
        N = NVAL[IN];

        for (IMAT = 1; IMAT <= NTYPES; IMAT++) {
          // 30

          // Do the tests only if DOTYPE[ IMAT ] is true.

          if (!DOTYPE[IMAT]) continue;

          // Test DGGRQF

          // Set up parameters with DLATB9 and generate test
          // matrices A and B with DLATMS.

          dlatb9(
            'GRQ',
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

          // Generate P by N matrix B

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

          NT = 4;

          dgrqts(
            M,
            P,
            N,
            A.asMatrix(LDA),
            AF.asMatrix(LDA),
            AQ.asMatrix(LDA),
            AR.asMatrix(LDA),
            LDA,
            TAUA,
            B.asMatrix(LDB),
            BF.asMatrix(LDB),
            BZ.asMatrix(LDB),
            BT.asMatrix(LDB),
            BWK.asMatrix(LDB),
            LDB,
            TAUB,
            WORK,
            LWORK,
            RWORK,
            RESULT,
          );

          // Print information about the tests that did not
          // pass the threshold.

          for (I = 1; I <= NT; I++) {
            // 10
            if (RESULT[I] >= THRESH) {
              if (NFAIL == 0 && FIRSTT) {
                FIRSTT = false;
                alahdg(NOUT, 'GRQ');
              }
              NOUT.println(
                ' M=${M.i4} P=${P.i4}, N=${N.i4}, type ${IMAT.i2}, test ${I.i2}, ratio=${RESULT[I].g13_6}',
              );
              NFAIL = NFAIL + 1;
            }
          } // 10
          NRUN = NRUN + NT;

          // Test DGGQRF

          // Set up parameters with DLATB9 and generate test
          // matrices A and B with DLATMS.

          dlatb9(
            'GQR',
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

          // Generate N-by-M matrix  A

          dlatms(
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

          // Generate N-by-P matrix  B

          dlatms(
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

          NT = 4;

          dgqrts(
            N,
            M,
            P,
            A,
            AF,
            AQ,
            AR,
            LDA,
            TAUA,
            B,
            BF,
            BZ,
            BT,
            BWK,
            LDB,
            TAUB,
            WORK,
            LWORK,
            RWORK,
            RESULT,
          );

          // Print information about the tests that did not
          // pass the threshold.

          for (I = 1; I <= NT; I++) {
            // 20
            if (RESULT[I] >= THRESH) {
              if (NFAIL == 0 && FIRSTT) {
                FIRSTT = false;
                alahdg(NOUT, PATH);
              }
              NOUT.println(
                ' N=${N.i4} M=${M.i4}, P=${P.i4}, type ${IMAT.i2}, test ${I.i2}, ratio=${RESULT[I].g13_6}',
              );
              NFAIL = NFAIL + 1;
            }
          } // 20
          NRUN = NRUN + NT;
        } // 30
      } // 40
    } // 50
  } // 60

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, 0);
}

void print9999(final Nout NOUT, final int info) {
  NOUT.println(' DLATMS in DCKGQR:    INFO = ${info.i5}');
}
