import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dsysv_aa.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../matgen/dlatms.dart';
import 'aladhd.dart';
import 'alaerh.dart';
import 'alasvm.dart';
import 'common.dart';
import 'derrvxx.dart';
import 'dlarhs.dart';
import 'dlatb4.dart';
import 'dpot02.dart';
import 'dsyt01_aa.dart';
import 'xlaenv.dart';

void ddrvsy_aa(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NRHS,
  final double THRESH,
  final bool TSTERR,
  final int NMAX,
  final Array<double> A_,
  final Array<double> AFAC_,
  final Array<double> AINV_,
  final Array<double> B_,
  final Array<double> X_,
  final Array<double> XACT_,
  final Array<double> WORK_,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Nout NOUT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DOTYPE = DOTYPE_.having();
  final NVAL = NVAL_.having();
  final A = A_.having();
  final AFAC = AFAC_.having();
  final AINV = AINV_.having();
  final B = B_.having();
  final X = X_.having();
  final XACT = XACT_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0;
  const NTYPES = 10, NTESTS = 3;
  const NFACT = 2;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L']; // FACTS = ['F', 'N'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  // Test path

  final PATH = '${'Double precision'[0]}SA';

  // Path to generate matrices

  final MATPATH = '${'Double precision'[0]}SY';

  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    // 10
    ISEED[I] = ISEEDY[I];
  } // 10

  // Test the error exits

  if (TSTERR) derrvx(PATH, NOUT);
  infoc.INFOT = 0;

  // Set the block size and minimum block size for testing.

  const NB = 1;
  const NBMIN = 2;
  xlaenv(1, NB);
  xlaenv(2, NBMIN);

  // Do for each value of N in NVAL

  for (var IN = 1; IN <= NN; IN++) {
    // 180
    final N = NVAL[IN];
    final LWORK = max(max(3 * N - 2, N * (1 + NB)), 1);
    final LDA = max(N, 1);
    var XTYPE = 'N';
    final NIMAT = N <= 0 ? 1 : NTYPES;

    for (var IMAT = 1; IMAT <= NIMAT; IMAT++) {
      // 170

      // Do the tests only if DOTYPE( IMAT ) is true.

      if (!DOTYPE[IMAT]) continue;

      // Skip types 3, 4, 5, or 6 if the matrix size is too small.

      final ZEROT = IMAT >= 3 && IMAT <= 6;
      if (ZEROT && N < IMAT - 2) continue;

      // Do first for UPLO = 'U', then for UPLO = 'L'

      for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
        // 160
        final UPLO = UPLOS[IUPLO - 1];

        // Set up parameters with DLATB4 and generate a test matrix
        // with DLATMS.

        final (:TYPE, :KL, :KU, :ANORM, :MODE, COND: CNDNUM, :DIST) =
            dlatb4(MATPATH, IMAT, N, N);

        srnamc.SRNAMT = 'DLATMS';
        dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
            UPLO, A.asMatrix(), LDA, WORK, INFO);

        // Check error code from DLATMS.

        if (INFO.value != 0) {
          alaerh(PATH, 'DLATMS', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
              NFAIL, NERRS, NOUT);
          continue;
        }

        // For types 3-6, zero one or more rows and columns of the
        // matrix to test that INFO.value is returned correctly.

        int IZERO;
        if (ZEROT) {
          if (IMAT == 3) {
            IZERO = 1;
          } else if (IMAT == 4) {
            IZERO = N;
          } else {
            IZERO = N ~/ 2 + 1;
          }

          if (IMAT < 6) {
            // Set row and column IZERO to zero.

            if (IUPLO == 1) {
              var IOFF = (IZERO - 1) * LDA;
              for (var I = 1; I <= IZERO - 1; I++) {
                // 20
                A[IOFF + I] = ZERO;
              } // 20
              IOFF = IOFF + IZERO;
              for (var I = IZERO; I <= N; I++) {
                // 30
                A[IOFF] = ZERO;
                IOFF = IOFF + LDA;
              } // 30
            } else {
              var IOFF = IZERO;
              for (var I = 1; I <= IZERO - 1; I++) {
                // 40
                A[IOFF] = ZERO;
                IOFF = IOFF + LDA;
              } // 40
              IOFF = IOFF - IZERO;
              for (var I = IZERO; I <= N; I++) {
                // 50
                A[IOFF + I] = ZERO;
              } // 50
            }
          } else {
            var IOFF = 0;
            if (IUPLO == 1) {
              // Set the first IZERO rows and columns to zero.

              for (var J = 1; J <= N; J++) {
                // 70
                final I2 = min(J, IZERO);
                for (var I = 1; I <= I2; I++) {
                  // 60
                  A[IOFF + I] = ZERO;
                } // 60
                IOFF = IOFF + LDA;
              } // 70
              IZERO = 1;
            } else {
              // Set the last IZERO rows and columns to zero.

              for (var J = 1; J <= N; J++) {
                // 90
                final I1 = max(J, IZERO);
                for (var I = I1; I <= N; I++) {
                  // 80
                  A[IOFF + I] = ZERO;
                } // 80
                IOFF = IOFF + LDA;
              } // 90
            }
          }
        } else {
          IZERO = 0;
        }

        for (var IFACT = 1; IFACT <= NFACT; IFACT++) {
          // 150

          // Do first for FACT = 'F', then for other values.

          // final FACT = FACTS[IFACT - 1];

          // Form an exact solution and set the right hand side.

          srnamc.SRNAMT = 'DLARHS';
          dlarhs(MATPATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A.asMatrix(),
              LDA, XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
          XTYPE = 'C';

          // --- Test DSYSV_AA  ---

          if (IFACT == 2) {
            dlacpy(UPLO, N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);
            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

            // Factor the matrix and solve the system using DSYSV_AA.

            srnamc.SRNAMT = 'DSYSV_AA';
            dsysv_aa(UPLO, N, NRHS, AFAC.asMatrix(), LDA, IWORK, X.asMatrix(),
                LDA, WORK, LWORK, INFO);

            // Adjust the expected value of INFO to account for
            // pivoting.

            var K = 0;
            if (IZERO > 0) {
              var J = 1;
              K = IZERO;
              while (true) {
                if (J == K) {
                  K = IWORK[J];
                } else if (IWORK[J] == K) {
                  K = J;
                }
                if (J < K) {
                  J++;
                  continue;
                }
                break;
              }
            }

            // Check error code from DSYSV_AA .

            if (INFO.value != K) {
              alaerh(PATH, 'DSYSV_AA ', INFO.value, K, UPLO, N, N, -1, -1, NRHS,
                  IMAT, NFAIL, NERRS, NOUT);
            } else if (INFO.value != 0) {
              //
            } else {
              // Reconstruct matrix from factors and compute
              // residual.

              dsyt01_aa(UPLO, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA, IWORK,
                  AINV.asMatrix(), LDA, RWORK, RESULT(1));

              // Compute residual of the computed solution.

              dlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
              dpot02(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                  WORK.asMatrix(), LDA, RWORK, RESULT(2));
              const NT = 2;

              // Print information about the tests that did not pass
              // the threshold.

              for (var K = 1; K <= NT; K++) {
                // 110
                if (RESULT[K] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                  NOUT.println(
                      ' DSYSV_AA , UPLO=\'${UPLO.a1}\', N =${N.i5}, type ${IMAT.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}');
                  NFAIL++;
                }
              } // 110
              NRUN += NT;
            } // 120
          }
        } // 150
      } // 160
    } // 170
  } // 180

  // Print a summary of the results.

  alasvm(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
