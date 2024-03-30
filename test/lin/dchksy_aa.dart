import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dsytrf_aa.dart';
import 'package:lapack/src/dsytrs_aa.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../matgen/dlatms.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'derrsyx.dart';
import 'dlarhs.dart';
import 'dlatb4.dart';
import 'dpot02.dart';
import 'dsyt01_aa.dart';
import 'xlaenv.dart';

void dchksy_aa(
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

  const ZERO = 0.0;
  const NTYPES = 10;
  const NTESTS = 9;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'];
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
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) derrsy(PATH, NOUT);
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
    final XTYPE = 'N';
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

        // Set up parameters with DLATB4 for the matrix generator
        // based on the type of matrix to be generated.

        final (:TYPE, :KL, :KU, :ANORM, :MODE, COND: CNDNUM, :DIST) =
            dlatb4(MATPATH, IMAT, N, N);

        // Generate a matrix with DLATMS.

        srnamc.SRNAMT = 'DLATMS';
        dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
            UPLO, A.asMatrix(), LDA, WORK, INFO);

        // Check error code from DLATMS and handle error.

        if (INFO.value != 0) {
          alaerh(PATH, 'DLATMS', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
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
                A[IOFF + I] = ZERO;
              }
              IOFF += IZERO;
              for (var I = IZERO; I <= N; I++) {
                A[IOFF] = ZERO;
                IOFF += LDA;
              }
            } else {
              var IOFF = IZERO;
              for (var I = 1; I <= IZERO - 1; I++) {
                A[IOFF] = ZERO;
                IOFF += LDA;
              }
              IOFF -= IZERO;
              for (var I = IZERO; I <= N; I++) {
                A[IOFF + I] = ZERO;
              }
            }
          } else {
            if (IUPLO == 1) {
              // Set the first IZERO rows and columns to zero.

              var IOFF = 0;
              for (var J = 1; J <= N; J++) {
                final I2 = min(J, IZERO);
                for (var I = 1; I <= I2; I++) {
                  A[IOFF + I] = ZERO;
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
                  A[IOFF + I] = ZERO;
                }
                IOFF += LDA;
              }
            }
          }
        } else {
          IZERO = 0;
        }

        // End generate the test matrix A.

        // Do for each value of NB in NBVAL

        for (var INB = 1; INB <= NNB; INB++) {
          // Set the optimal blocksize, which will be later
          // returned by ILAENV.

          final NB = NBVAL[INB];
          xlaenv(1, NB);

          // Copy the test matrix A into matrix AFAC which
          // will be factorized in place. This is needed to
          // preserve the test matrix A for subsequent tests.

          dlacpy(UPLO, N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);

          // Compute the L*D*L**T or U*D*U**T factorization of the
          // matrix. IWORK stores details of the interchanges and
          // the block structure of D. AINV is a work array for
          // block factorization, LWORK is the length of AINV.

          srnamc.SRNAMT = 'DSYTRF_AA';
          var LWORK = max(1, N * NB + N);
          dsytrf_aa(UPLO, N, AFAC.asMatrix(), LDA, IWORK, AINV, LWORK, INFO);

          // Adjust the expected value of INFO.value to account for
          // pivoting.

          // IF( IZERO > 0 ) THEN
          //    J = 1
          //    K = IZERO
          // 100 CONTINUE
          //    IF( J == K ) THEN
          //       K = IWORK( J )
          //    ELSE IF( IWORK( J ) == K ) THEN
          //       K = J
          //    END IF
          //    IF( J < K ) THEN
          //       J++
          //       GO TO 100
          //    END IF
          // ELSE
          var K = 0;
          // END IF

          // Check error code from DSYTRF and handle error.

          if (INFO.value != K) {
            alaerh(PATH, 'DSYTRF_AA', INFO.value, K, UPLO, N, N, -1, -1, NB,
                IMAT, NFAIL, NERRS, NOUT);
          }

          // +    TEST 1
          // Reconstruct matrix from factors and compute residual.

          dsyt01_aa(UPLO, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA, IWORK,
              AINV.asMatrix(), LDA, RWORK, RESULT(1));
          const NT = 1;

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

              srnamc.SRNAMT = 'DLARHS';
              dlarhs(
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
              dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

              srnamc.SRNAMT = 'DSYTRS_AA';
              LWORK = max(1, 3 * N - 2);
              dsytrs_aa(UPLO, N, NRHS, AFAC.asMatrix(), LDA, IWORK,
                  X.asMatrix(), LDA, WORK, LWORK, INFO);

              // Check error code from DSYTRS and handle error.

              if (INFO.value != 0) {
                if (IZERO == 0) {
                  alaerh(PATH, 'DSYTRS_AA', INFO.value, 0, UPLO, N, N, -1, -1,
                      NRHS, IMAT, NFAIL, NERRS, NOUT);
                }
              } else {
                dlacpy(
                    'Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);

                // Compute the residual for the solution

                dpot02(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
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
