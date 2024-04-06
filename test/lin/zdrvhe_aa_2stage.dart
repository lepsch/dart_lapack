import 'package:lapack/lapack.dart';

import '../matgen/zlatms.dart';
import 'aladhd.dart';
import 'alaerh.dart';
import 'alasvm.dart';
import 'common.dart';
import 'xlaenv.dart';
import 'zerrvx.dart';
import 'zlarhs.dart';
import 'zlatb4.dart';
import 'zpot02.dart';

void zdrvhe_aa_2stage(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NRHS,
  final double THRESH,
  final bool TSTERR,
  final int NMAX,
  final Array<Complex> A_,
  final Array<Complex> AFAC_,
  final Array<Complex> AINV_,
  final Array<Complex> B_,
  final Array<Complex> X_,
  final Array<Complex> XACT_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Nout NOUT,
) {
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

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NTYPES = 10, NTESTS = 3, NFACT = 2;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L']; // FACTS = ['F', 'N'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  // Test path

  final PATH = '${'Zomplex precision'[0]}H2';

  // Path to generate matrices

  final MATPATH = '${'Zomplex precision'[0]}HE';

  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) zerrvx(PATH, NOUT);
  infoc.INFOT = 0;

  // Set the block size and minimum block size for testing.

  final NB = 1;
  final NBMIN = 2;
  xlaenv(1, NB);
  xlaenv(2, NBMIN);

  // Do for each value of N in NVAL

  for (var IN = 1; IN <= NN; IN++) {
    final N = NVAL[IN];
    final LDA = max(N, 1);
    var XTYPE = 'N';
    final NIMAT = N <= 0 ? 1 : NTYPES;

    for (var IMAT = 1; IMAT <= NIMAT; IMAT++) {
      // Do the tests only if DOTYPE( IMAT ) is true.

      if (!DOTYPE[IMAT]) continue;

      // Skip types 3, 4, 5, or 6 if the matrix size is too small.

      final ZEROT = IMAT >= 3 && IMAT <= 6;
      if (ZEROT && N < IMAT - 2) continue;

      // Do first for UPLO = 'U', then for UPLO = 'L'

      for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
        final UPLO = UPLOS[IUPLO - 1];

        // Begin generate the test matrix A.

        // Set up parameters with ZLATB4 for the matrix generator
        // based on the type of matrix to be generated.

        final (:TYPE, :KL, :KU, :ANORM, :MODE, :CNDNUM, :DIST) =
            zlatb4(MATPATH, IMAT, N, N);

        // Generate a matrix with ZLATMS.

        srnamc.SRNAMT = 'ZLATMS';
        zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
            UPLO, A.asMatrix(), LDA, WORK, INFO);

        // Check error code from ZLATMS and handle error.

        if (INFO.value != 0) {
          alaerh(PATH, 'ZLATMS', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
              NFAIL, NERRS, NOUT);
          continue;
        }

        // For types 3-6, zero one or more rows and columns of
        // the matrix to test that INFO is returned correctly.

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
                A[IOFF + I] = Complex.zero;
              }
              IOFF += IZERO;
              for (var I = IZERO; I <= N; I++) {
                A[IOFF] = Complex.zero;
                IOFF += LDA;
              }
            } else {
              var IOFF = IZERO;
              for (var I = 1; I <= IZERO - 1; I++) {
                A[IOFF] = Complex.zero;
                IOFF += LDA;
              }
              IOFF -= IZERO;
              for (var I = IZERO; I <= N; I++) {
                A[IOFF + I] = Complex.zero;
              }
            }
          } else {
            var IOFF = 0;
            if (IUPLO == 1) {
              // Set the first IZERO rows and columns to zero.

              for (var J = 1; J <= N; J++) {
                final I2 = min(J, IZERO);
                for (var I = 1; I <= I2; I++) {
                  A[IOFF + I] = Complex.zero;
                }
                IOFF += LDA;
              }
              IZERO = 1;
            } else {
              // Set the first IZERO rows and columns to zero.

              var IOFF = 0;
              for (var J = 1; J <= N; J++) {
                final I1 = max(J, IZERO);
                for (var I = I1; I <= N; I++) {
                  A[IOFF + I] = Complex.zero;
                }
                IOFF += LDA;
              }
            }
          }
        } else {
          IZERO = 0;
        }

        // End generate the test matrix A.

        for (var IFACT = 1; IFACT <= NFACT; IFACT++) {
          // Do first for FACT = 'F', then for other values.

          // final FACT = FACTS[IFACT - 1];

          // Form an exact solution and set the right hand side.

          srnamc.SRNAMT = 'ZLARHS';
          zlarhs(MATPATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A.asMatrix(),
              LDA, XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
          XTYPE = 'C';

          // --- Test ZHESV_AA_2STAGE  ---

          if (IFACT == 2) {
            zlacpy(UPLO, N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);
            zlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

            // Factor the matrix and solve the system using ZHESV_AA.

            srnamc.SRNAMT = 'ZHESV_AA_2STAGE';
            final LWORK = min(max(1, N * NB), 3 * NMAX * NMAX);
            zhesv_aa_2stage(
                UPLO,
                N,
                NRHS,
                AFAC.asMatrix(),
                LDA,
                AINV,
                max(1, (3 * NB + 1) * N),
                IWORK,
                IWORK(1 + N),
                X.asMatrix(),
                LDA,
                WORK,
                LWORK,
                INFO);

            // Adjust the expected value of INFO to account for
            // pivoting.

            int K;
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
            } else {
              K = 0;
            }

            // Check error code from ZHESV_AA .

            if (INFO.value != K) {
              alaerh(PATH, 'ZHESV_AA', INFO.value, K, UPLO, N, N, -1, -1, NRHS,
                  IMAT, NFAIL, NERRS, NOUT);
            } else if (INFO.value != 0) {
              //
            } else {
              // Compute residual of the computed solution.

              zlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
              zpot02(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                  WORK.asMatrix(), LDA, RWORK, RESULT(1));

              // Reconstruct matrix from factors and compute
              // residual.

              // NEED TO CREATE ZHET01_AA_2STAGE
              //  CALL ZHET01_AA( UPLO, N, A, LDA, AFAC, LDA,
              //      IWORK, AINV, LDA, RWORK,
              //      RESULT( 2 ) )
              //  NT = 2
              const NT = 1;

              // Print information about the tests that did not pass
              // the threshold.

              for (var K = 1; K <= NT; K++) {
                if (RESULT[K] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                  NOUT.println(
                      ' ZHESV_AA_2STAGE, UPLO=\'${UPLO.a1}\', N =${N.i5}, type ${IMAT.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}');
                  NFAIL++;
                }
              }
              NRUN += NT;
            }
          }
        }
      }
    }
  }

  // Print a summary of the results.

  alasvm(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
