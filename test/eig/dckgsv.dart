import 'package:lapack/src/box.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/range.dart';
import 'package:test/test.dart';

import '../lin/alasum.dart';
import '../matgen/dlatms.dart';
import '../test_driver.dart';
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
  final TestDriver test,
  final String group,
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
  final DOTYPE = Array<bool>(NTYPES);
  final RESULT = Array<double>(NTESTS);
  const PATH = 'GSV';

  // Initialize constants and the random number seed.

  INFO.value = 0;
  var NRUN = 0;
  var NFAIL = 0;
  await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);
  final LDA = NMAX;
  final LDB = NMAX;
  final LDU = NMAX;
  final LDV = NMAX;
  final LDQ = NMAX;
  final LDR = NMAX;
  final LWORK = NMAX * NMAX;

  test.group(group, () {
    var FIRSTT = true;

    // Do for each value of M in MVAL.

    for (final IM in 1.through(NM)) {
      final M = MVAL[IM];
      final P = PVAL[IM];
      final N = NVAL[IM];

      for (final IMAT in 1.through(NTYPES)) {
        // Do the tests only if DOTYPE[ IMAT ] is true.
        final skip = !DOTYPE[IMAT];
        test('DCKGSV (M = $M, N = $N, P = $P TYPE = $IMAT)', () {
          final IINFO = Box(0);

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

          dlatms(M, N, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, ANORM, KLA,
              KUA, 'No packing', A.asMatrix(LDA), LDA, WORK, IINFO);
          test.expect(IINFO.value, 0);
          if (IINFO.value != 0) {
            print9999(NOUT, IINFO.value);
            INFO.value = IINFO.value.abs();
            return;
          }

          dlatms(P, N, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB, BNORM, KLB,
              KUB, 'No packing', B.asMatrix(LDB), LDB, WORK, IINFO);
          test.expect(IINFO.value, 0);
          if (IINFO.value != 0) {
            print9999(NOUT, IINFO.value);
            INFO.value = IINFO.value.abs();
            return;
          }

          const NT = 6;

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

          for (var I = 1; I <= NT; I++) {
            final reason =
                ' M=${M.i4} P=${P.i4}, N=${N.i4}, type ${IMAT.i2}, test ${I.i2}, ratio=${RESULT[I].g13_6}';
            test.expect(RESULT[I], lessThan(THRESH), reason: reason);
            if (RESULT[I] >= THRESH) {
              if (NFAIL == 0 && FIRSTT) {
                FIRSTT = false;
                alahdg(NOUT, PATH);
              }
              NOUT.println(reason);
              NFAIL++;
            }
          }
          NRUN += NT;
        }, skip: skip);
      }
    }
  });

  // Print a summary of the results.
  alasum(PATH, NOUT, NFAIL, NRUN, 0);
}

void print9999(final Nout NOUT, final int info) {
  NOUT.println(' DLATMS in DCKGSV   INFO = ${info.i5}');
}
