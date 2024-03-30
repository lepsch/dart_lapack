import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zhetrf_aa_2stage.dart';
import 'package:lapack/src/zhetrs_aa_2stage.dart';
import 'package:lapack/src/zlacpy.dart';

import '../matgen/zlatms.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'xlaenv.dart';
import 'zerrhe.dart';
import 'zlaipd.dart';
import 'zlarhs.dart';
import 'zlatb4.dart';
import 'zpot02.dart';

void zchkhe_aa_2stage(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NNB,
  final Array<int> NBVAL_,
  final int NNS,
  final Array<int> NSVAL_,
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
  final NBVAL = NBVAL_.having();
  final NSVAL = NSVAL_.having();
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
  const NTYPES = 10, NTESTS = 9;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'];
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

  if (TSTERR) zerrhe(PATH, NOUT);
  infoc.INFOT = 0;

  // Set the minimum block size for which the block routine should
  // be used, which will be later returned by ILAENV

  xlaenv(2, 2);

  // Do for each value of N in NVAL

  for (var IN = 1; IN <= NN; IN++) {
    final N = NVAL[IN];
    if (N > NMAX) {
      NFAIL++;
      NOUT.println(
          ' Invalid input value: ${'M '.a4}=${N.i6}; must be <=${NMAX.i6}');
      continue;
    }
    final LDA = max(N, 1);
    var XTYPE = 'N';
    final NIMAT = N <= 0 ? 1 : NTYPES;

    // Do for each value of matrix type IMAT

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

          // Skip all tests for this generated matrix

          continue;
        }

        // For matrix types 3-6, zero one or more rows and
        // columns of the matrix to test that INFO.value is returned
        // correctly.

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
            if (IUPLO == 1) {
              // Set the first IZERO rows and columns to zero.

              var IOFF = 0;
              for (var J = 1; J <= N; J++) {
                final I2 = min(J, IZERO);
                for (var I = 1; I <= I2; I++) {
                  A[IOFF + I] = Complex.zero;
                }
                IOFF += LDA;
              }
              IZERO = 1;
            } else {
              // Set the last IZERO rows and columns to zero.

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

        // End generate test matrix A.

        // Set the imaginary part of the diagonals.

        zlaipd(N, A, LDA + 1, 0);

        // Do for each value of NB in NBVAL

        for (var INB = 1; INB <= NNB; INB++) {
          // Set the optimal blocksize, which will be later
          // returned by ILAENV.

          final NB = NBVAL[INB];
          xlaenv(1, NB);

          // Copy the test matrix A into matrix AFAC which
          // will be factorized in place. This is needed to
          // preserve the test matrix A for subsequent tests.

          zlacpy(UPLO, N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);

          // Compute the L*D*L**T or U*D*U**T factorization of the
          // matrix. IWORK stores details of the interchanges and
          // the block structure of D. AINV is a work array for
          // block factorization, LWORK is the length of AINV.

          srnamc.SRNAMT = 'ZHETRF_AA_2STAGE';
          var LWORK = min(max(1, N * NB), 3 * NMAX * NMAX);
          zhetrf_aa_2stage(UPLO, N, AFAC.asMatrix(), LDA, AINV,
              max(1, (3 * NB + 1) * N), IWORK, IWORK(1 + N), WORK, LWORK, INFO);

          // Adjust the expected value of INFO.value to account for
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

          // Check error code from CHETRF and handle error.

          if (INFO.value != K) {
            alaerh(PATH, 'ZHETRF_AA_2STAGE', INFO.value, K, UPLO, N, N, -1, -1,
                NB, IMAT, NFAIL, NERRS, NOUT);
          }

          // +    TEST 1
          // Reconstruct matrix from factors and compute residual.

          // NEED TO CREATE ZHET01_AA_2STAGE
          //  CALL ZHET01_AA( UPLO, N, A, LDA, AFAC, LDA, IWORK,
          //      AINV, LDA, RWORK, RESULT( 1 ) )
          // NT = 1
          const NT = 0;

          // Print information about the tests that did not pass
          // the threshold.

          for (var K = 1; K <= NT; K++) {
            if (RESULT[K] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.println(
                  ' UPLO = \'${UPLO.a1}\', N =${N.i5}, NB =${NB.i4}, type ${IMAT.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}');
              NFAIL++;
            }
          }
          NRUN += NT;

          // Skip solver test if INFO.value is not 0.

          if (INFO.value == 0) {
            // Do for each value of NRHS in NSVAL.

            for (var IRHS = 1; IRHS <= NNS; IRHS++) {
              final NRHS = NSVAL[IRHS];

              // +    TEST 2 (Using TRS)
              // Solve and compute residual for  A * X = B.

              // Choose a set of NRHS random solution vectors
              // stored in XACT and set up the right hand side B

              srnamc.SRNAMT = 'ZLARHS';
              zlarhs(
                  MATPATH,
                  XTYPE,
                  UPLO,
                  ' ',
                  N,
                  N,
                  KL,
                  KU,
                  NRHS,
                  A.asMatrix(),
                  LDA,
                  XACT.asMatrix(),
                  LDA,
                  B.asMatrix(),
                  LDA,
                  ISEED,
                  INFO);
              zlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

              srnamc.SRNAMT = 'ZHETRS_AA_2STAGE';
              LWORK = max(1, 3 * N - 2);
              zhetrs_aa_2stage(
                  UPLO,
                  N,
                  NRHS,
                  AFAC.asMatrix(),
                  LDA,
                  AINV,
                  (3 * NB + 1) * N,
                  IWORK,
                  IWORK(1 + N),
                  X.asMatrix(),
                  LDA,
                  INFO);

              // Check error code from ZHETRS and handle error.

              if (INFO.value != 0) {
                if (IZERO == 0) {
                  alaerh(PATH, 'ZHETRS_AA_2STAGE', INFO.value, 0, UPLO, N, N,
                      -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT);
                }
              } else {
                zlacpy(
                    'Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);

                // Compute the residual for the solution

                zpot02(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                    WORK.asMatrix(), LDA, RWORK, RESULT(2));

                // Print information about the tests that did not pass
                // the threshold.

                for (var K = 2; K <= 2; K++) {
                  if (RESULT[K] >= THRESH) {
                    if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                    NOUT.println(
                        ' UPLO = \'${UPLO.a1}\', N =${N.i5}, NRHS=${NRHS.i3}, type ${IMAT.i2}, test(${K.i2}) =${RESULT[K].g12_5}');
                    NFAIL++;
                  }
                }
              }
              NRUN++;

              // End do for each value of NRHS in NSVAL.
            }
          }
        }
      }
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
