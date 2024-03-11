import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dsysv.dart';
import 'package:lapack/src/dsysvx.dart';
import 'package:lapack/src/dsysvxx.dart';
import 'package:lapack/src/dsytrf.dart';
import 'package:lapack/src/dsytri2.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../matgen/dlatms.dart';
import 'aladhd.dart';
import 'alaerh.dart';
import 'alasvm.dart';
import 'common.dart';
import 'debchvxx.dart';
import 'derrvxx.dart';
import 'dget04.dart';
import 'dget06.dart';
import 'dlarhs.dart';
import 'dlatb4.dart';
import 'dpot02.dart';
import 'dpot05.dart';
import 'dsyt01.dart';
import 'xlaenv.dart';

void ddrvsy(
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
  const ONE = 1.0, ZERO = 0.0;
  const NTYPES = 10, NTESTS = 6;
  const NFACT = 2;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS),
      BERR = Array<double>(NRHS),
      ERRBNDS_N = Matrix<double>(NRHS, 3),
      ERRBNDS_C = Matrix<double>(NRHS, 3);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'], FACTS = ['F', 'N'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}SY';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }
  var LWORK = max(2 * NMAX, NMAX * NRHS);

  // Test the error exits

  if (TSTERR) derrvx(PATH, NOUT);
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
        final UPLO = UPLOS[IUPLO];

        // Set up parameters with DLATB4 and generate a test matrix
        // with DLATMS.

        final (:TYPE, :KL, :KU, :ANORM, :MODE, COND: CNDNUM, :DIST) =
            dlatb4(PATH, IMAT, N, N);

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

        final int IZERO;
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
              IOFF = IOFF + IZERO;
              for (var I = IZERO; I <= N; I++) {
                A[IOFF] = ZERO;
                IOFF = IOFF + LDA;
              }
            } else {
              var IOFF = IZERO;
              for (var I = 1; I <= IZERO - 1; I++) {
                A[IOFF] = ZERO;
                IOFF = IOFF + LDA;
              }
              IOFF = IOFF - IZERO;
              for (var I = IZERO; I <= N; I++) {
                A[IOFF + I] = ZERO;
              }
            }
          } else {
            var IOFF = 0;
            if (IUPLO == 1) {
              // Set the first IZERO rows and columns to zero.

              for (var J = 1; J <= N; J++) {
                final I2 = min(J, IZERO);
                for (var I = 1; I <= I2; I++) {
                  A[IOFF + I] = ZERO;
                }
                IOFF = IOFF + LDA;
              }
            } else {
              // Set the last IZERO rows and columns to zero.

              for (var J = 1; J <= N; J++) {
                final I1 = max(J, IZERO);
                for (var I = I1; I <= N; I++) {
                  A[IOFF + I] = ZERO;
                }
                IOFF = IOFF + LDA;
              }
            }
          }
        } else {
          IZERO = 0;
        }

        var RCONDC = ZERO;
        for (var IFACT = 1; IFACT <= NFACT; IFACT++) {
          // Do first for FACT = 'F', then for other values.

          final FACT = FACTS[IFACT - 1];

          // Compute the condition number for comparison with
          // the value returned by DSYSVX.

          if (ZEROT) {
            if (IFACT == 1) continue;
            RCONDC = ZERO;
          } else if (IFACT == 1) {
            // Compute the 1-norm of A.

            final ANORM = dlansy('1', UPLO, N, A.asMatrix(), LDA, RWORK);

            // Factor the matrix A.

            dlacpy(UPLO, N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);
            dsytrf(UPLO, N, AFAC.asMatrix(), LDA, IWORK, WORK, LWORK, INFO);

            // Compute inv(A) and take its norm.

            dlacpy(UPLO, N, N, AFAC.asMatrix(), LDA, AINV.asMatrix(), LDA);
            LWORK = (N + NB + 1) * (NB + 3);
            dsytri2(UPLO, N, AINV.asMatrix(), LDA, IWORK, WORK, LWORK, INFO);
            final AINVNM = dlansy('1', UPLO, N, AINV.asMatrix(), LDA, RWORK);

            // Compute the 1-norm condition number of A.

            if (ANORM <= ZERO || AINVNM <= ZERO) {
              RCONDC = ONE;
            } else {
              RCONDC = (ONE / ANORM) / AINVNM;
            }
          }

          // Form an exact solution and set the right hand side.

          srnamc.SRNAMT = 'DLARHS';
          dlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A.asMatrix(), LDA,
              XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
          XTYPE = 'C';

          // --- Test DSYSV  ---

          if (IFACT == 2) {
            dlacpy(UPLO, N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);
            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

            // Factor the matrix and solve the system using DSYSV.

            srnamc.SRNAMT = 'DSYSV ';
            dsysv(UPLO, N, NRHS, AFAC.asMatrix(), LDA, IWORK, X.asMatrix(), LDA,
                WORK, LWORK, INFO);

            // Adjust the expected value of INFO.value to account for
            // pivoting.

            var K = IZERO;
            if (K > 0) {
              while (true) {
                if (IWORK[K] < 0) {
                  if (IWORK[K] != -K) {
                    K = -IWORK[K];
                    continue;
                  }
                } else if (IWORK[K] != K) {
                  K = IWORK[K];
                  continue;
                }
                break;
              }
            }

            // Check error code from DSYSV .

            if (INFO.value != K) {
              alaerh(PATH, 'DSYSV ', INFO.value, K, UPLO, N, N, -1, -1, NRHS,
                  IMAT, NFAIL, NERRS, NOUT);
            } else if (INFO.value != 0) {
              //
            } else {
              // Reconstruct matrix from factors and compute
              // residual.

              dsyt01(UPLO, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA, IWORK,
                  AINV.asMatrix(), LDA, RWORK, RESULT(1));

              // Compute residual of the computed solution.

              dlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
              dpot02(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                  WORK.asMatrix(), LDA, RWORK, RESULT(2));

              // Check solution from generated exact solution.

              dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                  RESULT(3));
              const NT = 3;

              // Print information about the tests that did not pass
              // the threshold.

              for (var K = 1; K <= NT; K++) {
                if (RESULT[K] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                  NOUT.println(
                      ' DSYSV , UPLO=\'${UPLO.a1}\', N =${N.i5}, type ${IMAT.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}');
                  NFAIL++;
                }
              }
              NRUN += NT;
            }
          }

          // --- Test DSYSVX ---

          if (IFACT == 2) dlaset(UPLO, N, N, ZERO, ZERO, AFAC.asMatrix(), LDA);
          dlaset('Full', N, NRHS, ZERO, ZERO, X.asMatrix(), LDA);

          // Solve the system and compute the condition number and
          // error bounds using DSYSVX.

          final RCOND = Box(0.0);
          srnamc.SRNAMT = 'DSYSVX';
          dsysvx(
              FACT,
              UPLO,
              N,
              NRHS,
              A.asMatrix(),
              LDA,
              AFAC.asMatrix(),
              LDA,
              IWORK,
              B.asMatrix(),
              LDA,
              X.asMatrix(),
              LDA,
              RCOND,
              RWORK,
              RWORK(NRHS + 1),
              WORK,
              LWORK,
              IWORK(N + 1),
              INFO);

          // Adjust the expected value of INFO.value to account for
          // pivoting.

          var K = IZERO;
          if (K > 0) {
            while (true) {
              if (IWORK[K] < 0) {
                if (IWORK[K] != -K) {
                  K = -IWORK[K];
                  continue;
                }
              } else if (IWORK[K] != K) {
                K = IWORK[K];
                continue;
              }
              break;
            }
          }

          // Check the error code from DSYSVX.

          if (INFO.value != K) {
            alaerh(PATH, 'DSYSVX', INFO.value, K, FACT + UPLO, N, N, -1, -1,
                NRHS, IMAT, NFAIL, NERRS, NOUT);
            continue;
          }

          final int K1;
          if (INFO.value == 0) {
            if (IFACT >= 2) {
              // Reconstruct matrix from factors and compute
              // residual.

              dsyt01(UPLO, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA, IWORK,
                  AINV.asMatrix(), LDA, RWORK(2 * NRHS + 1), RESULT(1));
              K1 = 1;
            } else {
              K1 = 2;
            }

            // Compute residual of the computed solution.

            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
            dpot02(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                WORK.asMatrix(), LDA, RWORK(2 * NRHS + 1), RESULT(2));

            // Check solution from generated exact solution.

            dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                RESULT(3));

            // Check the error bounds from iterative refinement.

            dpot05(
                UPLO,
                N,
                NRHS,
                A.asMatrix(),
                LDA,
                B.asMatrix(),
                LDA,
                X.asMatrix(),
                LDA,
                XACT.asMatrix(),
                LDA,
                RWORK,
                RWORK(NRHS + 1),
                RESULT(4));
          } else {
            K1 = 6;
          }

          // Compare RCOND from DSYSVX with the computed value
          // in RCONDC.

          RESULT[6] = dget06(RCOND.value, RCONDC);

          // Print information about the tests that did not pass
          // the threshold.

          for (var K = K1; K <= 6; K++) {
            if (RESULT[K] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
              NOUT.println(
                  ' DSYSVX, FACT=\'${FACT.a1}\', UPLO=\'${UPLO.a1}\', N =${N.i5}, type ${IMAT.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}');
              NFAIL++;
            }
          }
          NRUN += 7 - K1;

          // --- Test DSYSVXX ---

          // Restore the matrices A and B.

          if (IFACT == 2) dlaset(UPLO, N, N, ZERO, ZERO, AFAC.asMatrix(), LDA);
          dlaset('Full', N, NRHS, ZERO, ZERO, X.asMatrix(), LDA);

          // Solve the system and compute the condition number
          // and error bounds using DSYSVXX.

          srnamc.SRNAMT = 'DSYSVXX';
          const N_ERR_BNDS = 3;
          final EQUED = Box('N'), RPVGRW_SVXX = Box(0.0);
          dsysvxx(
              FACT,
              UPLO,
              N,
              NRHS,
              A.asMatrix(),
              LDA,
              AFAC.asMatrix(),
              LDA,
              IWORK,
              EQUED,
              WORK(N + 1),
              B.asMatrix(),
              LDA,
              X.asMatrix(),
              LDA,
              RCOND,
              RPVGRW_SVXX,
              BERR,
              N_ERR_BNDS,
              ERRBNDS_N,
              ERRBNDS_C,
              0,
              Array<double>(1),
              WORK,
              IWORK(N + 1),
              INFO);

          // Adjust the expected value of INFO.value to account for
          // pivoting.

          var K2 = IZERO;
          if (K2 > 0) {
            while (true) {
              if (IWORK[K2] < 0) {
                if (IWORK[K2] != -K2) {
                  K2 = -IWORK[K2];
                  continue;
                }
              } else if (IWORK[K2] != K2) {
                K2 = IWORK[K2];
                continue;
              }
              break;
            }
          }

          // Check the error code from DSYSVXX.

          if (INFO.value != K2 && INFO.value <= N) {
            alaerh(PATH, 'DSYSVXX', INFO.value, K2, FACT + UPLO, N, N, -1, -1,
                NRHS, IMAT, NFAIL, NERRS, NOUT);
            continue;
          }

          final int K3;
          if (INFO.value == 0) {
            if (IFACT >= 2) {
              // Reconstruct matrix from factors and compute
              // residual.

              dsyt01(UPLO, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA, IWORK,
                  AINV.asMatrix(), LDA, RWORK(2 * NRHS + 1), RESULT(1));
              K3 = 1;
            } else {
              K3 = 2;
            }

            // Compute residual of the computed solution.

            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
            dpot02(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                WORK.asMatrix(), LDA, RWORK(2 * NRHS + 1), RESULT(2));

            // Check solution from generated exact solution.

            dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                RESULT(3));

            // Check the error bounds from iterative refinement.

            dpot05(
                UPLO,
                N,
                NRHS,
                A.asMatrix(),
                LDA,
                B.asMatrix(),
                LDA,
                X.asMatrix(),
                LDA,
                XACT.asMatrix(),
                LDA,
                RWORK,
                RWORK(NRHS + 1),
                RESULT(4));
          } else {
            K3 = 6;
          }

          // Compare RCOND.value from DSYSVXX with the computed value
          // in RCONDC.

          RESULT[6] = dget06(RCOND.value, RCONDC);

          // Print information about the tests that did not pass
          // the threshold.

          for (var K = K3; K <= 6; K++) {
            if (RESULT[K] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
              NOUT.println(
                  ' DSYSVXX, FACT=\'${FACT.a1}\', UPLO=\'${UPLO.a1}\', N =${N.i5}, type ${IMAT.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}');
              NFAIL++;
            }
          }
          NRUN += 7 - K3;
        }
      }
    }
  }

  // Print a summary of the results.

  alasvm(PATH, NOUT, NFAIL, NRUN, NERRS.value);

  // Test Error Bounds from DSYSVXX

  debchvxx(THRESH, PATH);
}
