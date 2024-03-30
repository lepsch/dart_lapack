import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgels.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../matgen/dlatms.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'derrqr.dart';
import 'dgennd.dart';
import 'dget02.dart';
import 'dlarhs.dart';
import 'dlatb4.dart';
import 'dqrt01.dart';
import 'dqrt01p.dart';
import 'dqrt02.dart';
import 'dqrt03.dart';
import 'xlaenv.dart';

void dchkqr(
  final Array<bool> DOTYPE_,
  final int NM,
  final Array<int> MVAL_,
  final int NN,
  final Array<int> NVAL_,
  final int NNB,
  final Array<int> NBVAL_,
  final Array<int> NXVAL_,
  final int NRHS,
  final double THRESH,
  final bool TSTERR,
  final int NMAX,
  final Array<double> A_,
  final Array<double> AF_,
  final Array<double> AQ_,
  final Array<double> AR_,
  final Array<double> AC_,
  final Array<double> B_,
  final Array<double> X_,
  final Array<double> XACT_,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Nout NOUT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DOTYPE = DOTYPE_.having();
  final MVAL = MVAL_.having();
  final NVAL = NVAL_.having();
  final NBVAL = NBVAL_.having();
  final NXVAL = NXVAL_.having();
  final A = A_.having();
  final AF = AF_.having();
  final AQ = AQ_.having();
  final AR = AR_.having();
  final AC = AC_.having();
  final B = B_.having();
  final X = X_.having();
  final XACT = XACT_.having();
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  // final IWORK = IWORK_.having();
  const NTESTS = 9;
  const NTYPES = 8;
  const ZERO = 0.0;
  final ISEED = Array<int>(4), KVAL = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}QR';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) derrqr(PATH, NOUT);
  infoc.INFOT = 0;
  xlaenv(2, 2);

  final LDA = NMAX;
  final LWORK = NMAX * max(NMAX, NRHS).toInt();

  // Do for each value of M in MVAL.

  for (var IM = 1; IM <= NM; IM++) {
    final M = MVAL[IM];

    // Do for each value of N in NVAL.

    for (var IN = 1; IN <= NN; IN++) {
      final N = NVAL[IN];
      final MINMN = min(M, N);
      for (var IMAT = 1; IMAT <= NTYPES; IMAT++) {
        // Do the tests only if DOTYPE( IMAT ) is true.

        if (!DOTYPE[IMAT]) continue;

        // Set up parameters with DLATB4 and generate a test matrix
        // with DLATMS.

        final (:TYPE, :KL, :KU, :ANORM, :MODE, COND: CNDNUM, :DIST) =
            dlatb4(PATH, IMAT, M, N);

        srnamc.SRNAMT = 'DLATMS';
        dlatms(M, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
            'No packing', A.asMatrix(), LDA, WORK, INFO);

        // Check error code from DLATMS.

        if (INFO.value != 0) {
          alaerh(PATH, 'DLATMS', INFO.value, 0, ' ', M, N, -1, -1, -1, IMAT,
              NFAIL, NERRS, NOUT);
          continue;
        }

        // Set some values for K: the first value must be MINMN,
        // corresponding to the call of DQRT01; other values are
        // used in the calls of DQRT02, and must not exceed MINMN.

        KVAL[1] = MINMN;
        KVAL[2] = 0;
        KVAL[3] = 1;
        KVAL[4] = MINMN ~/ 2;
        final int NK;
        if (MINMN == 0) {
          NK = 1;
        } else if (MINMN == 1) {
          NK = 2;
        } else if (MINMN <= 3) {
          NK = 3;
        } else {
          NK = 4;
        }

        // Do for each value of K in KVAL

        for (var IK = 1; IK <= NK; IK++) {
          final K = KVAL[IK];

          // Do for each pair of values (NB,NX) in NBVAL and NXVAL.

          for (var INB = 1; INB <= NNB; INB++) {
            final NB = NBVAL[INB];
            xlaenv(1, NB);
            final NX = NXVAL[INB];
            xlaenv(3, NX);
            for (var I = 1; I <= NTESTS; I++) {
              RESULT[I] = ZERO;
            }
            // ignore: unused_local_variable
            var NT = 2;
            if (IK == 1) {
              // Test DGEQRF

              dqrt01(M, N, A.asMatrix(), AF.asMatrix(), AQ.asMatrix(),
                  AR.asMatrix(), LDA, TAU, WORK, LWORK, RWORK, RESULT(1));

              // Test DGEQRFP

              dqrt01p(M, N, A.asMatrix(), AF.asMatrix(), AQ.asMatrix(),
                  AR.asMatrix(), LDA, TAU, WORK, LWORK, RWORK, RESULT(8));
              if (!dgennd(M, N, AF.asMatrix(), LDA)) RESULT[9] = 2 * THRESH;
              NT++;
            } else if (M >= N) {
              // Test DORGQR, using factorization
              // returned by DQRT01

              dqrt02(M, N, K, A.asMatrix(), AF.asMatrix(), AQ.asMatrix(),
                  AR.asMatrix(), LDA, TAU, WORK, LWORK, RWORK, RESULT(1));
            }
            if (M >= K) {
              // Test DORMQR, using factorization returned
              // by DQRT01

              dqrt03(M, N, K, AF.asMatrix(), AC.asMatrix(), AR.asMatrix(),
                  AQ.asMatrix(), LDA, TAU, WORK, LWORK, RWORK, RESULT(3));
              NT += 4;

              // If M>=N and K=N, call DGELS to solve a system
              // with NRHS right hand sides and compute the
              // residual.

              if (K == N && INB == 1) {
                // Generate a solution and set the right
                // hand side.

                srnamc.SRNAMT = 'DLARHS';
                dlarhs(
                    PATH,
                    'New',
                    'Full',
                    'No transpose',
                    M,
                    N,
                    0,
                    0,
                    NRHS,
                    A.asMatrix(),
                    LDA,
                    XACT.asMatrix(),
                    LDA,
                    B.asMatrix(),
                    LDA,
                    ISEED,
                    INFO);

                dlacpy('Full', M, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

                // Reset AF. DGELS overwrites the matrix with
                // its factorization.

                dlacpy('Full', M, N, A.asMatrix(), LDA, AF.asMatrix(), LDA);

                srnamc.SRNAMT = 'DGELS';
                dgels('No transpose', M, N, NRHS, AF.asMatrix(), LDA,
                    X.asMatrix(), LDA, WORK, LWORK, INFO);

                // Check error code from DGELS.

                if (INFO.value != 0) {
                  alaerh(PATH, 'DGELS', INFO.value, 0, 'N', M, N, NRHS, -1, NB,
                      IMAT, NFAIL, NERRS, NOUT);
                }

                dget02('No transpose', M, N, NRHS, A.asMatrix(), LDA,
                    X.asMatrix(), LDA, B.asMatrix(), LDA, RWORK, RESULT(7));
                NT++;
              }
            }

            // Print information about the tests that did not
            // pass the threshold.

            for (var I = 1; I <= NTESTS; I++) {
              if (RESULT[I] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                NOUT.println(
                    ' M=${M.i5}, N=${N.i5}, K=${K.i5}, NB=${NB.i4}, NX=${NX.i5}, type ${IMAT.i2}, test(${I.i2})=${RESULT[I].g12_5}');
                NFAIL++;
              }
            }
            NRUN += NTESTS;
          }
        }
      }
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
