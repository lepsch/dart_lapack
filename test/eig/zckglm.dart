import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../matgen/zlarnd.dart';
import '../matgen/zlatms.dart';
import 'alahdg.dart';
import 'alareq.dart';
import 'alasum.dart';
import 'dlatb9.dart';
import 'zglmts.dart';

Future<void> zckglm(
  final int NN,
  final Array<int> NVAL_,
  final Array<int> MVAL_,
  final Array<int> PVAL_,
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
  final NVAL = NVAL_.having();
  final MVAL = MVAL_.having();
  final PVAL = PVAL_.having();
  final A = A_.having();
  final AF = AF_.having();
  final B = B_.having();
  final BF = BF_.having();
  final X = X_.having();
  final ISEED = ISEED_.having(length: 4);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();

  const NTYPES = 8;
  bool FIRSTT;
  String PATH;
  int I, IK, IMAT, LDA, LDB, LWORK, M, N, NFAIL, NRUN, P;

  final DOTYPE = Array<bool>(NTYPES);
  final DISTA = Box(''), DISTB = Box(''), TYPE = Box('');
  final ANORM = Box(0.0),
      BNORM = Box(0.0),
      CNDNMA = Box(0.0),
      CNDNMB = Box(0.0),
      RESID = Box(0.0);
  final IINFO = Box(0),
      KLA = Box(0),
      KLB = Box(0),
      KUA = Box(0),
      KUB = Box(0),
      MODEA = Box(0),
      MODEB = Box(0);

  // Initialize constants.

  PATH = 'GLM';
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
    // 10
    M = MVAL[IK];
    P = PVAL[IK];
    N = NVAL[IK];
    if (M > N || N > M + P) {
      if (FIRSTT) {
        NOUT.println();
        FIRSTT = false;
      }
      NOUT.println(
          ' *** Invalid input  for GLM:  M = ${M.i6}, P = ${P.i6}, N = ${N.i6};\n     must satisfy M <= N <= M+P  (this set of values will be skipped)');
    }
  } // 10
  FIRSTT = true;

  // Do for each value of M in MVAL.

  for (IK = 1; IK <= NN; IK++) {
    // 40
    M = MVAL[IK];
    P = PVAL[IK];
    N = NVAL[IK];
    if (M > N || N > M + P) continue;

    for (IMAT = 1; IMAT <= NTYPES; IMAT++) {
      // 30

      // Do the tests only if DOTYPE( IMAT ) is true.

      if (!DOTYPE[IMAT]) continue;

      // Set up parameters with DLATB9 and generate test
      // matrices A and B with ZLATMS.

      dlatb9(PATH, IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, ANORM, BNORM, MODEA,
          MODEB, CNDNMA, CNDNMB, DISTA, DISTB);

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

      // Generate random left hand side vector of GLM

      for (I = 1; I <= N; I++) {
        // 20
        X[I] = zlarnd(2, ISEED);
      } // 20

      zglmts(
          N,
          M,
          P,
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
          WORK,
          LWORK,
          RWORK,
          RESID);

      // Print information about the tests that did not
      // pass the threshold.

      if (RESID.value >= THRESH) {
        if (NFAIL == 0 && FIRSTT) {
          FIRSTT = false;
          alahdg(NOUT, PATH);
        }
        NOUT.println(
            ' N=${N.i4} M=${M.i4}, P=${P.i4}, type ${IMAT.i2}, test ${I.i2}, ratio=${RESID.value.g13_6}');
        NFAIL++;
      }
      NRUN++;
    } // 30
  } // 40

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, 0);

//  9999 FORMAT( '${.i5}');
//  9998 FORMAT( ;
//  9997 FORMAT(  )
}

void _print9999(Nout nout, int info) {
  nout.println(' ZLATMS in ZCKGLM INFO = ${info.i5}');
}
