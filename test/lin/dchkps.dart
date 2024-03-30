import 'dart:math';

import 'package:lapack/lapack.dart';
import 'package:test/test.dart';

import '../matgen/dlatmt.dart';
import '../test_driver.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'derrps.dart';
import 'dlatb5.dart';
import 'dpst01.dart';
import 'xlaenv.dart';

void dchkps(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NNB,
  final Array<int> NBVAL_,
  final int NRANK,
  final Array<int> RANKVAL_,
  final double THRESH,
  final bool TSTERR,
  final int NMAX,
  final Array<double> A_,
  final Array<double> AFAC_,
  final Array<double> PERM_,
  final Array<int> PIV_,
  final Array<double> WORK_,
  final Array<double> RWORK_,
  final Nout NOUT,
  final TestDriver test,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DOTYPE = DOTYPE_.having();
  final NVAL = NVAL_.having();
  final NBVAL = NBVAL_.having();
  final RANKVAL = RANKVAL_.having();
  final PIV = PIV_.having();
  final A = A_.having();
  final AFAC = AFAC_.having();
  final PERM = PERM_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();

  const ONE = 1.0;
  const NTYPES = 9;
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'];

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}PS';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  final ISEED = Array.fromList(ISEEDY);

  test.group('error exits', () {
    // Test the error exits
    if (TSTERR) derrps(PATH, NOUT, test);
    test.tearDown(() {
      infoc.INFOT = 0;
    });
  });

  test.setUp(() {
    xlaenv(2, 2);
  });

  // Do for each value of N in NVAL

  for (final IN in 1.through(NN)) {
    final N = NVAL[IN];
    final LDA = max(N, 1);
    final NIMAT = N <= 0 ? 1 : NTYPES;

    for (final IMAT in 1.through(NIMAT)) {
      // Do the tests only if DOTYPE( IMAT ) is true.
      final skip = !DOTYPE[IMAT];

      test('DCHKPS (IN=$IN IMAT=$IMAT)', () {
        final INFO = Box(0);

        // Do for each value of RANK in RANKVAL
        for (var IRANK = 1; IRANK <= NRANK; IRANK++) {
          // Only repeat test 3 to 5 for different ranks
          // Other tests use full rank
          if ((IMAT < 3 || IMAT > 5) && IRANK > 1) continue;

          final RANK = ((N * RANKVAL[IRANK]) / 100.0).ceil();

          // Do first for UPLO = 'U', then for UPLO = 'L'

          for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
            final UPLO = UPLOS[IUPLO - 1];

            // Set up parameters with DLATB5 and generate a test matrix
            // with DLATMT.
            final (:TYPE, :KL, :KU, :ANORM, :MODE, :CNDNUM, :DIST) =
                dlatb5(PATH, IMAT, N);

            srnamc.SRNAMT = 'DLATMT';
            dlatmt(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, RANK,
                KL, KU, UPLO, A.asMatrix(), LDA, WORK, INFO);

            // Check error code from DLATMT.
            test.expect(INFO.value, 0);
            if (INFO.value != 0) {
              alaerh(PATH, 'DLATMT', INFO.value, 0, UPLO, N, N, -1, -1, -1,
                  IMAT, NFAIL, NERRS, NOUT);
              continue;
            }

            // Do for each value of NB in NBVAL

            for (var INB = 1; INB <= NNB; INB++) {
              final NB = NBVAL[INB];
              xlaenv(1, NB);

              // Compute the pivoted L*L' or U'*U factorization
              // of the matrix.

              dlacpy(UPLO, N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);
              srnamc.SRNAMT = 'DPSTRF';

              // Use default tolerance
              const TOL = -ONE;
              final COMPRANK = Box(0);
              dpstrf(UPLO, N, AFAC.asMatrix(), LDA, PIV, COMPRANK, TOL, WORK,
                  INFO);

              // Check error code from DPSTRF.
              const IZERO = 0;
              // `alaerh` ignores INFO == 0
              if (INFO.value != 0) {
                test.expect(INFO.value, isNot(lessThan(IZERO)));
                if (RANK == N) {
                  test.expect(INFO.value, isNot(IZERO));
                } else if (RANK < N) {
                  test.expect(INFO.value, isNot(lessThanOrEqualTo(IZERO)));
                }
              }
              if ((INFO.value < IZERO) ||
                  (INFO.value != IZERO && RANK == N) ||
                  (INFO.value <= IZERO && RANK < N)) {
                alaerh(PATH, 'DPSTRF', INFO.value, IZERO, UPLO, N, N, -1, -1,
                    NB, IMAT, NFAIL, NERRS, NOUT);
                continue;
              }

              // Skip the test if INFO is not 0.
              if (INFO.value != 0) continue;

              // Reconstruct matrix from factors and compute residual.

              // PERM holds permuted L*L^T or U^T*U

              final RESULT = Box(0.0);
              dpst01(UPLO, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA,
                  PERM.asMatrix(), LDA, PIV, RWORK, RESULT, COMPRANK.value);

              // Print information about the tests that did not pass
              // the threshold or where computed rank was not RANK.

              if (N == 0) COMPRANK.value = 0;
              final RANKDIFF = RANK - COMPRANK.value;
              final reason =
                  ' UPLO = \'${UPLO.a1}\', N =${N.i5}, RANK =${RANK.i3}, Diff =${RANKDIFF.i5}, NB =${NB.i4}, type ${IMAT.i2}, Ratio =${RESULT.value.g12_5}';
              test.expect(RESULT.value, lessThan(THRESH));
              if (RESULT.value >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                NOUT.println(reason);
                NFAIL++;
              }
              NRUN++;
            }
          }
        }
      }, skip: skip);
    }
  }

  // Print a summary of the results.
  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
