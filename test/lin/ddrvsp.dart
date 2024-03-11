import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlansp.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dspsv.dart';
import 'package:lapack/src/dspsvx.dart';
import 'package:lapack/src/dsptrf.dart';
import 'package:lapack/src/dsptri.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../matgen/dlatms.dart';
import 'aladhd.dart';
import 'alaerh.dart';
import 'alasvm.dart';
import 'common.dart';
import 'derrvxx.dart';
import 'dget04.dart';
import 'dget06.dart';
import 'dlarhs.dart';
import 'dlatb4.dart';
import 'dppt02.dart';
import 'dppt05.dart';
import 'dspt01.dart';

void ddrvsp(
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
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const FACTS = ['F', 'N'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}SP';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I];
  }

  // Test the error exits

  if (TSTERR) derrvx(PATH, NOUT);
  infoc.INFOT = 0;

  // Do for each value of N in NVAL

  for (var IN = 1; IN <= NN; IN++) {
    final N = NVAL[IN];
    final LDA = max(N, 1);
    final NPP = N * (N + 1) ~/ 2;
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
        final (UPLO, PACKIT) = IUPLO == 1 ? ('U', 'C') : ('L', 'R');

        // Set up parameters with DLATB4 and generate a test matrix
        // with DLATMS.

        final (:TYPE, :KL, :KU, :ANORM, :MODE, COND: CNDNUM, :DIST) =
            dlatb4(PATH, IMAT, N, N);

        srnamc.SRNAMT = 'DLATMS';
        dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
            PACKIT, A.asMatrix(), LDA, WORK, INFO);

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
              var IOFF = (IZERO - 1) * IZERO ~/ 2;
              for (var I = 1; I <= IZERO - 1; I++) {
                A[IOFF + I] = ZERO;
              }
              IOFF = IOFF + IZERO;
              for (var I = IZERO; I <= N; I++) {
                A[IOFF] = ZERO;
                IOFF = IOFF + I;
              }
            } else {
              var IOFF = IZERO;
              for (var I = 1; I <= IZERO - 1; I++) {
                A[IOFF] = ZERO;
                IOFF = IOFF + N - I;
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
                IOFF = IOFF + J;
              }
            } else {
              // Set the last IZERO rows and columns to zero.

              for (var J = 1; J <= N; J++) {
                final I1 = max(J, IZERO);
                for (var I = I1; I <= N; I++) {
                  A[IOFF + I] = ZERO;
                }
                IOFF = IOFF + N - J;
              }
            }
          }
        } else {
          IZERO = 0;
        }

        var RCONDC = 0.0;
        for (var IFACT = 1; IFACT <= NFACT; IFACT++) {
          // Do first for FACT = 'F', then for other values.

          final FACT = FACTS[IFACT - 1];

          // Compute the condition number for comparison with
          // the value returned by DSPSVX.

          if (ZEROT) {
            if (IFACT == 1) continue;
            RCONDC = ZERO;
          } else if (IFACT == 1) {
            // Compute the 1-norm of A.

            final ANORM = dlansp('1', UPLO, N, A, RWORK);

            // Factor the matrix A.

            dcopy(NPP, A, 1, AFAC, 1);
            dsptrf(UPLO, N, AFAC, IWORK, INFO);

            // Compute inv(A) and take its norm.

            dcopy(NPP, AFAC, 1, AINV, 1);
            dsptri(UPLO, N, AINV, IWORK, WORK, INFO);
            final AINVNM = dlansp('1', UPLO, N, AINV, RWORK);

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

          // --- Test DSPSV  ---

          if (IFACT == 2) {
            dcopy(NPP, A, 1, AFAC, 1);
            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

            // Factor the matrix and solve the system using DSPSV.

            srnamc.SRNAMT = 'DSPSV ';
            dspsv(UPLO, N, NRHS, AFAC, IWORK, X.asMatrix(), LDA, INFO);

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

            // Check error code from DSPSV .

            if (INFO.value != K) {
              alaerh(PATH, 'DSPSV ', INFO.value, K, UPLO, N, N, -1, -1, NRHS,
                  IMAT, NFAIL, NERRS, NOUT);
            } else if (INFO.value != 0) {
              //
            } else {
              // Reconstruct matrix from factors and compute
              // residual.

              dspt01(UPLO, N, A, AFAC, IWORK, AINV.asMatrix(), LDA, RWORK,
                  RESULT(1));

              // Compute residual of the computed solution.

              dlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
              dppt02(UPLO, N, NRHS, A, X.asMatrix(), LDA, WORK.asMatrix(), LDA,
                  RWORK, RESULT(2));

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
                      ' DSPSV , UPLO=\'${UPLO.a1}\', N =${N.i5}, type ${IMAT.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}');
                  NFAIL++;
                }
              }
              NRUN += NT;
            }
          }

          // --- Test DSPSVX ---

          if (IFACT == 2 && NPP > 0) {
            dlaset('Full', NPP, 1, ZERO, ZERO, AFAC.asMatrix(), NPP);
          }
          dlaset('Full', N, NRHS, ZERO, ZERO, X.asMatrix(), LDA);

          // Solve the system and compute the condition number and
          // error bounds using DSPSVX.

          final RCOND = Box(0.0);
          srnamc.SRNAMT = 'DSPSVX';
          dspsvx(
              FACT,
              UPLO,
              N,
              NRHS,
              A,
              AFAC,
              IWORK,
              B.asMatrix(),
              LDA,
              X.asMatrix(),
              LDA,
              RCOND,
              RWORK,
              RWORK(NRHS + 1),
              WORK,
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

          // Check the error code from DSPSVX.

          if (INFO.value != K) {
            alaerh(PATH, 'DSPSVX', INFO.value, K, FACT + UPLO, N, N, -1, -1,
                NRHS, IMAT, NFAIL, NERRS, NOUT);
            continue;
          }

          final int K1;
          if (INFO.value == 0) {
            if (IFACT >= 2) {
              // Reconstruct matrix from factors and compute
              // residual.

              dspt01(UPLO, N, A, AFAC, IWORK, AINV.asMatrix(), LDA,
                  RWORK(2 * NRHS + 1), RESULT(1));
              K1 = 1;
            } else {
              K1 = 2;
            }

            // Compute residual of the computed solution.

            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
            dppt02(UPLO, N, NRHS, A, X.asMatrix(), LDA, WORK.asMatrix(), LDA,
                RWORK(2 * NRHS + 1), RESULT(2));

            // Check solution from generated exact solution.

            dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                RESULT(3));

            // Check the error bounds from iterative refinement.

            dppt05(UPLO, N, NRHS, A, B.asMatrix(), LDA, X.asMatrix(), LDA,
                XACT.asMatrix(), LDA, RWORK, RWORK(NRHS + 1), RESULT(4));
          } else {
            K1 = 6;
          }

          // Compare RCOND from DSPSVX with the computed value
          // in RCONDC.

          RESULT[6] = dget06(RCOND.value, RCONDC);

          // Print information about the tests that did not pass
          // the threshold.

          for (var K = K1; K <= 6; K++) {
            if (RESULT[K] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
              NOUT.println(
                  ' DSPSVX, FACT=\'${FACT.a1}\', UPLO=\'${UPLO.a1}\', N =${N.i5}, type ${IMAT.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}');
              NFAIL++;
            }
          }
          NRUN += 7 - K1;
        }
      }
    }
  }

  // Print a summary of the results.

  alasvm(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
