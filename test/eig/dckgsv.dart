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
  final Array<int> MVAL_,
  final Array<int> PVAL_,
  final Array<int> NVAL_,
  final int NMATS,
  final Array<int> ISEED_,
  final double THRESH,
  final int NMAX,
  final Array<double> A_,
  final Array<double> AF_,
  final Array<double> B_,
  final Array<double> BF_,
  final Array<double> U_,
  final Array<double> V_,
  final Array<double> Q_,
  final Array<double> ALPHA_,
  final Array<double> BETA_,
  final Array<double> R_,
  final Array<int> IWORK_,
  final Array<double> WORK_,
  final Array<double> RWORK_,
  final Nin NIN,
  final Nout NOUT,
  final Box<int> INFO,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final MVAL = MVAL_.having();
  final PVAL = PVAL_.having();
  final NVAL = NVAL_.having();
  final ISEED = ISEED_.having();
  final A = A_.having();
  final AF = AF_.having();
  final B = B_.having();
  final BF = BF_.having();
  final U = U_.having();
  final V = V_.having();
  final Q = Q_.having();
  final ALPHA = ALPHA_.having();
  final BETA = BETA_.having();
  final R = R_.having();
  final IWORK = IWORK_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
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
  final IINFO = Box(0);

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

      dlatms(M, N, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, ANORM, KLA, KUA,
          'No packing', A.asMatrix(LDA), LDA, WORK, IINFO);
      if (IINFO.value != 0) {
        print9999(NOUT, IINFO.value);
        INFO.value = IINFO.value.abs();
        continue;
      }

      dlatms(P, N, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB, BNORM, KLB, KUB,
          'No packing', B.asMatrix(LDB), LDB, WORK, IINFO);
      if (IINFO.value != 0) {
        print9999(NOUT, IINFO.value);
        INFO.value = IINFO.value.abs();
        continue;
      }

      NT = 6;

      dgsvts3(
          M,
          P,
          N,
          A.asMatrix(LDA),
          AF.asMatrix(LDA),
          LDA,
          B.asMatrix(LDB),
          BF.asMatrix(LDB),
          LDB,
          U.asMatrix(LDU),
          LDU,
          V.asMatrix(LDV),
          LDV,
          Q.asMatrix(LDQ),
          LDQ,
          ALPHA,
          BETA,
          R.asMatrix(LDR),
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

void print9999(final Nout NOUT, final int info) {
  NOUT.println(' DLATMS in DCKGSV   INFO = ${info.i5}');
}
