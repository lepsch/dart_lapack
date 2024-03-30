import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlansy.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zsysv.dart';
import 'package:lapack/src/zsysvx.dart';
import 'package:lapack/src/zsytrf.dart';
import 'package:lapack/src/zsytri2.dart';

import '../matgen/zlatms.dart';
import 'aladhd.dart';
import 'alaerh.dart';
import 'alasvm.dart';
import 'common.dart';
import 'dget06.dart';
import 'xlaenv.dart';
import 'zerrvxx.dart';
import 'zget04.dart';
import 'zlarhs.dart';
import 'zlatb4.dart';
import 'zlatsy.dart';
import 'zpot05.dart';
import 'zsyt01.dart';
import 'zsyt02.dart';

void zdrvsy(
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
  const NTYPES = 11, NTESTS = 6, NFACT = 2;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'], FACTS = ['F', 'N'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Zomplex precision'[0]}SY';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }
  var LWORK = max(2 * NMAX, NMAX * NRHS);

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

        int IZERO = 0, KL = 0, KU = 0;
        if (IMAT != NTYPES) {
          // Set up parameters with ZLATB4 and generate a test
          // matrix with ZLATMS.

          final (:TYPE, KL: KLTMP, KU: KUTMP, :ANORM, :MODE, :CNDNUM, :DIST) =
              zlatb4(PATH, IMAT, N, N);
          (KL, KU) = (KLTMP, KUTMP);

          srnamc.SRNAMT = 'ZLATMS';
          zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
              UPLO, A.asMatrix(), LDA, WORK, INFO);

          // Check error code from ZLATMS.

          if (INFO.value != 0) {
            alaerh(PATH, 'ZLATMS', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
                NFAIL, NERRS, NOUT);
            continue;
          }

          // For types 3-6, zero one or more rows and columns of
          // the matrix to test that INFO is returned correctly.

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
                // Set the first IZERO rows to zero.

                var IOFF = 0;
                for (var J = 1; J <= N; J++) {
                  final I2 = min(J, IZERO);
                  for (var I = 1; I <= I2; I++) {
                    A[IOFF + I] = Complex.zero;
                  }
                  IOFF += LDA;
                }
              } else {
                // Set the last IZERO rows to zero.

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
        } else {
          // IMAT = NTYPES:  Use a special block diagonal matrix to
          // test alternate code for the 2-by-2 blocks.

          zlatsy(UPLO, N, A.asMatrix(), LDA, ISEED);
        }

        var RCONDC = ZERO;
        for (var IFACT = 1; IFACT <= NFACT; IFACT++) {
          // Do first for FACT = 'F', then for other values.

          final FACT = FACTS[IFACT - 1];

          // Compute the condition number for comparison with
          // the value returned by ZSYSVX.

          if (ZEROT) {
            if (IFACT == 1) continue;
            RCONDC = ZERO;
          } else if (IFACT == 1) {
            // Compute the 1-norm of A.

            final ANORM = zlansy('1', UPLO, N, A.asMatrix(), LDA, RWORK);

            // Factor the matrix A.

            zlacpy(UPLO, N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);
            zsytrf(UPLO, N, AFAC.asMatrix(), LDA, IWORK, WORK, LWORK, INFO);

            // Compute inv(A) and take its norm.

            zlacpy(UPLO, N, N, AFAC.asMatrix(), LDA, AINV.asMatrix(), LDA);
            LWORK = (N + NB + 1) * (NB + 3);
            zsytri2(UPLO, N, AINV.asMatrix(), LDA, IWORK, WORK, LWORK, INFO);
            final AINVNM = zlansy('1', UPLO, N, AINV.asMatrix(), LDA, RWORK);

            // Compute the 1-norm condition number of A.

            if (ANORM <= ZERO || AINVNM <= ZERO) {
              RCONDC = ONE;
            } else {
              RCONDC = (ONE / ANORM) / AINVNM;
            }
          }

          // Form an exact solution and set the right hand side.

          srnamc.SRNAMT = 'ZLARHS';
          zlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A.asMatrix(), LDA,
              XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
          XTYPE = 'C';

          // --- Test ZSYSV  ---

          if (IFACT == 2) {
            zlacpy(UPLO, N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);
            zlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

            // Factor the matrix and solve the system using ZSYSV.

            srnamc.SRNAMT = 'ZSYSV ';
            zsysv(UPLO, N, NRHS, AFAC.asMatrix(), LDA, IWORK, X.asMatrix(), LDA,
                WORK, LWORK, INFO);

            // Adjust the expected value of INFO to account for
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

            // Check error code from ZSYSV .

            if (INFO.value != K) {
              alaerh(PATH, 'ZSYSV ', INFO.value, K, UPLO, N, N, -1, -1, NRHS,
                  IMAT, NFAIL, NERRS, NOUT);
            } else if (INFO.value != 0) {
              //
            } else {
              // Reconstruct matrix from factors and compute
              // residual.

              zsyt01(UPLO, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA, IWORK,
                  AINV.asMatrix(), LDA, RWORK, RESULT(1));

              // Compute residual of the computed solution.

              zlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
              zsyt02(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                  WORK.asMatrix(), LDA, RWORK, RESULT(2));

              // Check solution from generated exact solution.

              zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                  RESULT(3));
              const NT = 3;

              // Print information about the tests that did not pass
              // the threshold.

              for (var K = 1; K <= NT; K++) {
                if (RESULT[K] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                  NOUT.println(
                      ' ZSYSV , UPLO=\'${UPLO.a1}\', N =${N.i5}, type ${IMAT.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}');
                  NFAIL++;
                }
              }
              NRUN += NT;
            }
          }

          // --- Test ZSYSVX ---

          if (IFACT == 2) {
            zlaset(
                UPLO, N, N, Complex.zero, Complex.zero, AFAC.asMatrix(), LDA);
          }
          zlaset(
              'Full', N, NRHS, Complex.zero, Complex.zero, X.asMatrix(), LDA);

          // Solve the system and compute the condition number and
          // error bounds using ZSYSVX.

          final RCOND = Box(ZERO);
          srnamc.SRNAMT = 'ZSYSVX';
          zsysvx(
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
              RWORK(2 * NRHS + 1),
              INFO);

          // Adjust the expected value of INFO to account for
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

          // Check the error code from ZSYSVX.

          if (INFO.value != K) {
            alaerh(PATH, 'ZSYSVX', INFO.value, K, FACT + UPLO, N, N, -1, -1,
                NRHS, IMAT, NFAIL, NERRS, NOUT);
          } else {
            final int K1;
            if (INFO.value == 0) {
              if (IFACT >= 2) {
                // Reconstruct matrix from factors and compute
                // residual.

                zsyt01(UPLO, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA, IWORK,
                    AINV.asMatrix(), LDA, RWORK(2 * NRHS + 1), RESULT(1));
                K1 = 1;
              } else {
                K1 = 2;
              }

              // Compute residual of the computed solution.

              zlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
              zsyt02(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                  WORK.asMatrix(), LDA, RWORK(2 * NRHS + 1), RESULT(2));

              // Check solution from generated exact solution.

              zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                  RESULT(3));

              // Check the error bounds from iterative refinement.

              zpot05(
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

            // Compare RCOND from ZSYSVX with the computed value
            // in RCONDC.

            RESULT[6] = dget06(RCOND.value, RCONDC);

            // Print information about the tests that did not pass
            // the threshold.

            for (var K = K1; K <= 6; K++) {
              if (RESULT[K] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                NOUT.println(
                    ' ZSYSVX, FACT=\'${FACT.a1}\', UPLO=\'${UPLO.a1}\', N =${N.i5}, type ${IMAT.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}');
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
}
