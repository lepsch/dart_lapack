import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../lin/alasum.dart';
import '../matgen/dlatms.dart';
import 'alahdg.dart';
import 'alareq.dart';
import 'dlarhs.dart';
import 'dlatb9.dart';
import 'dlsets.dart';

Future<void> dcklse(
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
  final MVAL = MVAL_.having();
  final PVAL = PVAL_.having();
  final NVAL = NVAL_.having();
  final ISEED = ISEED_.having();
  final A = A_.having();
  final AF = AF_.having();
  final B = B_.having();
  final BF = BF_.having();
  final X = X_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const NTESTS = 7;
  const NTYPES = 8;
  bool FIRSTT;
  int I, IK, IMAT, LDA, LDB, LWORK, M, N, NFAIL, NRUN, NT, P;
  final DOTYPE = Array<bool>(NTYPES);
  final RESULT = Array<double>(NTESTS);
  final DISTA = Box(''), DISTB = Box(''), TYPE = Box('');
  final IINFO = Box(0),
      KLA = Box(0),
      KLB = Box(0),
      KUA = Box(0),
      KUB = Box(0),
      MODEA = Box(0),
      MODEB = Box(0);
  final ANORM = Box(0.0),
      BNORM = Box(0.0),
      CNDNMA = Box(0.0),
      CNDNMB = Box(0.0);
  const PATH = 'LSE';

  // Initialize constants and the random number seed.

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
      // matrices A and B with DLATMS.

      dlatb9(PATH, IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, ANORM, BNORM, MODEA,
          MODEB, CNDNMA, CNDNMB, DISTA, DISTB);

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
          IINFO);
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
          IINFO);
      if (IINFO.value != 0) {
        print9999(NOUT, IINFO.value);
        INFO.value = (IINFO.value).abs();
        continue;
      }

      // Generate the right-hand sides C and D for the LSE.
      var LD1 = max(N, 1), LD2 = max(M, 1);

      dlarhs(
          'DGE',
          'New solution',
          'Upper',
          'N',
          M,
          N,
          max(M - 1, 0),
          max(N - 1, 0),
          1,
          A.asMatrix(LDA),
          LDA,
          X(4 * NMAX + 1).asMatrix(LD1),
          LD1,
          X.asMatrix(LD2),
          LD2,
          ISEED,
          IINFO);

      LD1 = max(N, 1);
      LD2 = max(P, 1);
      dlarhs(
          'DGE',
          'Computed',
          'Upper',
          'N',
          P,
          N,
          max(P - 1, 0),
          max(N - 1, 0),
          1,
          B.asMatrix(LDB),
          LDB,
          X(4 * NMAX + 1).asMatrix(LD1),
          LD1,
          X(2 * NMAX + 1).asMatrix(LD2),
          LD2,
          ISEED,
          IINFO);

      NT = 2;

      dlsets(
          M,
          P,
          N,
          A.asMatrix(LDA),
          AF.asMatrix(LDA),
          LDA,
          B.asMatrix(LDB),
          BF.asMatrix(LDB),
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

void print9999(final Nout nout, final int info) {
  nout.println(' DLATMS in DCKLSE   INFO = ${info.i5}');
}
