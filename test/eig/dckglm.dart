import 'package:lapack/src/box.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../lin/alasum.dart';
import '../matgen/dlarnd.dart';
import '../matgen/dlatms.dart';
import 'alahdg.dart';
import 'alareq.dart';
import 'dglmts.dart';
import 'dlatb9.dart';

Future<void> dckglm(
  final int NN,
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
  final Array<double> X_,
  final Array<double> WORK_,
  final Array<double> RWORK_,
  final Nin NIN,
  final Nout NOUT,
  final Box<int> INFO,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final MVAL = MVAL_.dim();
  final PVAL = PVAL_.dim();
  final NVAL = NVAL_.dim();
  final ISEED = ISEED_.dim();
  final A = A_.dim();
  final AF = AF_.dim();
  final B = B_.dim();
  final BF = BF_.dim();
  final X = X_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
  const NTYPES = 8;
  bool FIRSTT;
  final TYPE = Box(''), DISTA = Box(''), DISTB = Box('');
  int I, IK, IMAT, LDA, LDB, LWORK, M, N, NFAIL, NRUN, P;
  final KLA = Box(0),
      KLB = Box(0),
      KUA = Box(0),
      KUB = Box(0),
      MODEA = Box(0),
      MODEB = Box(0);
  final ANORM = Box(0.0),
      BNORM = Box(0.0),
      CNDNMA = Box(0.0),
      CNDNMB = Box(0.0),
      RESID = Box(0.0);
  final DOTYPE = Array<bool>(NTYPES);
  final IINFO = Box(0);
  const PATH = 'GLM';

  // Initialize constants.

  INFO.value = 0;
  NRUN = 0;
  NFAIL = 0;
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
    if (M > N || N > M + P) {
      if (FIRSTT) {
        NOUT.println();
        FIRSTT = false;
      }
      NOUT.println(
          ' *** Invalid input  for GLM:  M = ${M.i6}, P = ${P.i6}, N = ${N.i6};\n'
          '     must satisfy M <= N <= M+P  (this set of values will be skipped)');
    }
  }
  FIRSTT = true;

  // Do for each value of M in MVAL.

  for (IK = 1; IK <= NN; IK++) {
    M = MVAL[IK];
    P = PVAL[IK];
    N = NVAL[IK];
    if (M > N || N > M + P) continue;

    for (IMAT = 1; IMAT <= NTYPES; IMAT++) {
      // Do the tests only if DOTYPE[ IMAT ] is true.
      if (!DOTYPE[IMAT]) continue;

      // Set up parameters with DLATB9 and generate test
      // matrices A and B with DLATMS.

      dlatb9(PATH, IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, ANORM, BNORM, MODEA,
          MODEB, CNDNMA, CNDNMB, DISTA, DISTB);

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
          IINFO);
      if (IINFO.value != 0) {
        print9999(NOUT, IINFO.value);
        INFO.value = (IINFO.value).abs();
        continue;
      }

      dlatms(
          N,
          P,
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
          IINFO);
      if (IINFO.value != 0) {
        print9999(NOUT, IINFO.value);
        INFO.value = (IINFO.value).abs();
        continue;
      }

      // Generate random left hand side vector of GLM

      for (I = 1; I <= N; I++) {
        X[I] = dlarnd(2, ISEED);
      }

      dglmts(N, M, P, A.asMatrix(LDA), AF.asMatrix(LDA), LDA, B.asMatrix(LDB), BF.asMatrix(LDB), LDB, X, X(NMAX + 1), X(2 * NMAX + 1),
          X(3 * NMAX + 1), WORK, LWORK, RWORK, RESID);

      // Print information about the tests that did not
      // pass the threshold.

      if (RESID.value >= THRESH) {
        if (NFAIL == 0 && FIRSTT) {
          FIRSTT = false;
          alahdg(NOUT, PATH);
        }
        NOUT.println(
            ' N=${M.i4} M=${N.i4}, P=${P.i4}, type ${IMAT.i2}, test 1, ratio=${RESID.value.g13_6}');
        NFAIL = NFAIL + 1;
      }
      NRUN = NRUN + 1;
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, 0);
}

void print9999(final Nout NOUT, final int info) {
  NOUT.println(' DLATMS in DCKGLM INFO = ${info.i5}');
}
