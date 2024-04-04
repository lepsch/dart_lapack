import 'dart:math';

import 'package:lapack/lapack.dart';

import '../matgen/zlatms.dart';
import 'aladhd.dart';
import 'alaerh.dart';
import 'alasvm.dart';
import 'common.dart';
import 'dget06.dart';
import 'xlaenv.dart';
import 'zebchvxx.dart';
import 'zerrvx.dart';
import 'zget04.dart';
import 'zlaipd.dart';
import 'zlarhs.dart';
import 'zlatb4.dart';
import 'zpot01.dart';
import 'zpot02.dart';
import 'zpot05.dart';

void zdrvpo(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NRHS,
  final double THRESH,
  final bool TSTERR,
  final int NMAX,
  final Array<Complex> A_,
  final Array<Complex> AFAC_,
  final Array<Complex> ASAV_,
  final Array<Complex> B_,
  final Array<Complex> BSAV_,
  final Array<Complex> X_,
  final Array<Complex> XACT_,
  final Array<double> S_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Nout NOUT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DOTYPE = DOTYPE_.having();
  final NVAL = NVAL_.having();
  final A = A_.having();
  final AFAC = AFAC_.having();
  final ASAV = ASAV_.having();
  final B = B_.having();
  final BSAV = BSAV_.having();
  final X = X_.having();
  final XACT = XACT_.having();
  final S = S_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  const NTYPES = 9, NTESTS = 6;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS),
      BERR = Array<double>(NRHS),
      ERRBNDS_N = Matrix<double>(NRHS, 3),
      ERRBNDS_C = Matrix<double>(NRHS, 3);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'];
  const FACTS = ['F', 'N', 'E'];
  const EQUEDS = ['N', 'Y'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Zomplex precision'[0]}PO';
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

      // Skip types 3, 4, or 5 if the matrix size is too small.

      final ZEROT = IMAT >= 3 && IMAT <= 5;
      if (ZEROT && N < IMAT - 2) continue;

      // Do first for UPLO = 'U', then for UPLO = 'L'

      for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
        final UPLO = UPLOS[IUPLO - 1];

        // Set up parameters with ZLATB4 and generate a test matrix
        // with ZLATMS.

        final (:TYPE, :KL, :KU, :ANORM, :MODE, :CNDNUM, :DIST) =
            zlatb4(PATH, IMAT, N, N);

        srnamc.SRNAMT = 'ZLATMS';
        zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
            UPLO, A.asMatrix(), LDA, WORK, INFO);

        // Check error code from ZLATMS.

        if (INFO.value != 0) {
          alaerh(PATH, 'ZLATMS', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
              NFAIL, NERRS, NOUT);
          continue;
        }

        // For types 3-5, zero one row and column of the matrix to
        // test that INFO is returned correctly.

        final int IZERO;
        if (ZEROT) {
          if (IMAT == 3) {
            IZERO = 1;
          } else if (IMAT == 4) {
            IZERO = N;
          } else {
            IZERO = N ~/ 2 + 1;
          }
          var IOFF = (IZERO - 1) * LDA;

          // Set row and column IZERO of A to 0.

          if (IUPLO == 1) {
            for (var I = 1; I <= IZERO - 1; I++) {
              A[IOFF + I] = Complex.zero;
            }
            IOFF += IZERO;
            for (var I = IZERO; I <= N; I++) {
              A[IOFF] = Complex.zero;
              IOFF += LDA;
            }
          } else {
            IOFF = IZERO;
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
          IZERO = 0;
        }

        // Set the imaginary part of the diagonals.

        zlaipd(N, A, LDA + 1, 0);

        // Save a copy of the matrix A in ASAV.

        zlacpy(UPLO, N, N, A.asMatrix(), LDA, ASAV.asMatrix(), LDA);

        for (var IEQUED = 1; IEQUED <= 2; IEQUED++) {
          final EQUED = Box(EQUEDS[IEQUED - 1]);
          final NFACT = IEQUED == 1 ? 3 : 1;

          var RCONDC = ZERO, ROLDC = ZERO;
          final SCOND = Box(ZERO), AMAX = Box(ZERO);
          for (var IFACT = 1; IFACT <= NFACT; IFACT++) {
            final FACT = FACTS[IFACT - 1];
            final PREFAC = lsame(FACT, 'F');
            final NOFACT = lsame(FACT, 'N');
            final EQUIL = lsame(FACT, 'E');

            if (ZEROT) {
              if (PREFAC) continue;
              RCONDC = ZERO;
            } else if (!lsame(FACT, 'N')) {
              // Compute the condition number for comparison with
              // the value returned by ZPOSVX (FACT = 'N' reuses
              // the condition number from the previous iteration
              // with FACT = 'F').

              zlacpy(UPLO, N, N, ASAV.asMatrix(), LDA, AFAC.asMatrix(), LDA);
              if (EQUIL || IEQUED > 1) {
                // Compute row and column scale factors to
                // equilibrate the matrix A.

                zpoequ(N, AFAC.asMatrix(), LDA, S, SCOND, AMAX, INFO);
                if (INFO.value == 0 && N > 0) {
                  if (IEQUED > 1) SCOND.value = ZERO;

                  // Equilibrate the matrix.

                  zlaqhe(UPLO, N, AFAC.asMatrix(), LDA, S, SCOND.value,
                      AMAX.value, EQUED);
                }
              }

              // Save the condition number of the
              // non-equilibrated system for use in ZGET04.

              if (EQUIL) ROLDC = RCONDC;

              // Compute the 1-norm of A.

              final ANORM = zlanhe('1', UPLO, N, AFAC.asMatrix(), LDA, RWORK);

              // Factor the matrix A.

              zpotrf(UPLO, N, AFAC.asMatrix(), LDA, INFO);

              // Form the inverse of A.

              zlacpy(UPLO, N, N, AFAC.asMatrix(), LDA, A.asMatrix(), LDA);
              zpotri(UPLO, N, A.asMatrix(), LDA, INFO);

              // Compute the 1-norm condition number of A.

              final AINVNM = zlanhe('1', UPLO, N, A.asMatrix(), LDA, RWORK);
              if (ANORM <= ZERO || AINVNM <= ZERO) {
                RCONDC = ONE;
              } else {
                RCONDC = (ONE / ANORM) / AINVNM;
              }
            }

            // Restore the matrix A.

            zlacpy(UPLO, N, N, ASAV.asMatrix(), LDA, A.asMatrix(), LDA);

            // Form an exact solution and set the right hand side.

            srnamc.SRNAMT = 'ZLARHS';
            zlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A.asMatrix(),
                LDA, XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
            XTYPE = 'C';
            zlacpy('Full', N, NRHS, B.asMatrix(), LDA, BSAV.asMatrix(), LDA);

            if (NOFACT) {
              // --- Test ZPOSV  ---

              // Compute the L*L' or U'*U factorization of the
              // matrix and solve the system.

              zlacpy(UPLO, N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);
              zlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

              srnamc.SRNAMT = 'ZPOSV';
              zposv(
                  UPLO, N, NRHS, AFAC.asMatrix(), LDA, X.asMatrix(), LDA, INFO);

              // Check error code from ZPOSV .

              if (INFO.value != IZERO) {
                alaerh(PATH, 'ZPOSV ', INFO.value, IZERO, UPLO, N, N, -1, -1,
                    NRHS, IMAT, NFAIL, NERRS, NOUT);
              } else if (INFO.value != 0) {
                //
              } else {
                // Reconstruct matrix from factors and compute
                // residual.

                zpot01(UPLO, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA, RWORK,
                    RESULT(1));

                // Compute residual of the computed solution.

                zlacpy(
                    'Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
                zpot02(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                    WORK.asMatrix(), LDA, RWORK, RESULT(2));

                // Check solution from generated exact solution.

                zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                    RESULT(3));
                const NT = 3;

                // Print information about the tests that did not
                // pass the threshold.

                for (var K = 1; K <= NT; K++) {
                  if (RESULT[K] >= THRESH) {
                    if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                    NOUT.println(
                        ' ZPOSV , UPLO=\'${UPLO.a1}\', N =${N.i5}, type ${IMAT.i1}, test(${K.i1})=${RESULT[K].g12_5}');
                    NFAIL++;
                  }
                }
                NRUN += NT;
              }
            }

            {
              // --- Test ZPOSVX ---

              if (!PREFAC) {
                zlaset(UPLO, N, N, Complex.zero, Complex.zero, AFAC.asMatrix(),
                    LDA);
              }
              zlaset('Full', N, NRHS, Complex.zero, Complex.zero, X.asMatrix(),
                  LDA);
              if (IEQUED > 1 && N > 0) {
                // Equilibrate the matrix if FACT='F' and
                // EQUED='Y'.

                zlaqhe(UPLO, N, A.asMatrix(), LDA, S, SCOND.value, AMAX.value,
                    EQUED);
              }

              // Solve the system and compute the condition number
              // and error bounds using ZPOSVX.

              srnamc.SRNAMT = 'ZPOSVX';
              final RCOND = Box(ZERO);
              zposvx(
                  FACT,
                  UPLO,
                  N,
                  NRHS,
                  A.asMatrix(),
                  LDA,
                  AFAC.asMatrix(),
                  LDA,
                  EQUED,
                  S,
                  B.asMatrix(),
                  LDA,
                  X.asMatrix(),
                  LDA,
                  RCOND,
                  RWORK,
                  RWORK(NRHS + 1),
                  WORK,
                  RWORK(2 * NRHS + 1),
                  INFO);

              // Check the error code from ZPOSVX.

              if (INFO.value != IZERO) {
                alaerh(PATH, 'ZPOSVX', INFO.value, IZERO, FACT + UPLO, N, N, -1,
                    -1, NRHS, IMAT, NFAIL, NERRS, NOUT);
                continue;
              }

              final int K1;
              if (INFO.value == 0) {
                if (!PREFAC) {
                  // Reconstruct matrix from factors and compute
                  // residual.

                  zpot01(UPLO, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA,
                      RWORK(2 * NRHS + 1), RESULT(1));
                  K1 = 1;
                } else {
                  K1 = 2;
                }

                // Compute residual of the computed solution.

                zlacpy('Full', N, NRHS, BSAV.asMatrix(), LDA, WORK.asMatrix(),
                    LDA);
                zpot02(UPLO, N, NRHS, ASAV.asMatrix(), LDA, X.asMatrix(), LDA,
                    WORK.asMatrix(), LDA, RWORK(2 * NRHS + 1), RESULT(2));

                // Check solution from generated exact solution.

                if (NOFACT || (PREFAC && lsame(EQUED.value, 'N'))) {
                  zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                      RCONDC, RESULT(3));
                } else {
                  zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                      ROLDC, RESULT(3));
                }

                // Check the error bounds from iterative
                // refinement.

                zpot05(
                    UPLO,
                    N,
                    NRHS,
                    ASAV.asMatrix(),
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

              // Compare RCOND from ZPOSVX with the computed value
              // in RCONDC.

              RESULT[6] = dget06(RCOND.value, RCONDC);

              // Print information about the tests that did not pass
              // the threshold.

              for (var K = K1; K <= 6; K++) {
                if (RESULT[K] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                  if (PREFAC) {
                    NOUT.print9997('ZPOSVX', FACT, UPLO, N, EQUED.value, IMAT,
                        K, RESULT[K]);
                  } else {
                    NOUT.print9998('ZPOSVX', FACT, UPLO, N, IMAT, K, RESULT[K]);
                  }
                  NFAIL++;
                }
              }
              NRUN += 7 - K1;
            }

            {
              // --- Test ZPOSVXX ---

              // Restore the matrices A and B.

              zlacpy('Full', N, N, ASAV.asMatrix(), LDA, A.asMatrix(), LDA);
              zlacpy('Full', N, NRHS, BSAV.asMatrix(), LDA, B.asMatrix(), LDA);
              if (!PREFAC) {
                zlaset(UPLO, N, N, Complex.zero, Complex.zero, AFAC.asMatrix(),
                    LDA);
              }
              zlaset('Full', N, NRHS, Complex.zero, Complex.zero, X.asMatrix(),
                  LDA);
              if (IEQUED > 1 && N > 0) {
                // Equilibrate the matrix if FACT='F' and
                // EQUED='Y'.

                zlaqhe(UPLO, N, A.asMatrix(), LDA, S, SCOND.value, AMAX.value,
                    EQUED);
              }

              // Solve the system and compute the condition number
              // and error bounds using ZPOSVXX.

              srnamc.SRNAMT = 'ZPOSVXX';
              const N_ERR_BNDS = 3;
              final RPVGRW_SVXX = Box(ZERO), RCOND = Box(ZERO);
              zposvxx(
                  FACT,
                  UPLO,
                  N,
                  NRHS,
                  A.asMatrix(),
                  LDA,
                  AFAC.asMatrix(),
                  LDA,
                  EQUED,
                  S,
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
                  RWORK(2 * NRHS + 1),
                  INFO);

              // Check the error code from ZPOSVXX.

              if (INFO.value == N + 1) continue;
              if (INFO.value != IZERO) {
                alaerh(PATH, 'ZPOSVXX', INFO.value, IZERO, FACT + UPLO, N, N,
                    -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT);
                continue;
              }

              final int K1;
              if (INFO.value == 0) {
                if (!PREFAC) {
                  // Reconstruct matrix from factors and compute
                  // residual.

                  zpot01(UPLO, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA,
                      RWORK(2 * NRHS + 1), RESULT(1));
                  K1 = 1;
                } else {
                  K1 = 2;
                }

                // Compute residual of the computed solution.

                zlacpy('Full', N, NRHS, BSAV.asMatrix(), LDA, WORK.asMatrix(),
                    LDA);
                zpot02(UPLO, N, NRHS, ASAV.asMatrix(), LDA, X.asMatrix(), LDA,
                    WORK.asMatrix(), LDA, RWORK(2 * NRHS + 1), RESULT(2));

                // Check solution from generated exact solution.

                if (NOFACT || (PREFAC && lsame(EQUED.value, 'N'))) {
                  zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                      RCONDC, RESULT(3));
                } else {
                  zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                      ROLDC, RESULT(3));
                }

                // Check the error bounds from iterative
                // refinement.

                zpot05(
                    UPLO,
                    N,
                    NRHS,
                    ASAV.asMatrix(),
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

              // Compare RCOND from ZPOSVXX with the computed value
              // in RCONDC.

              RESULT[6] = dget06(RCOND.value, RCONDC);

              // Print information about the tests that did not pass
              // the threshold.

              for (var K = K1; K <= 6; K++) {
                if (RESULT[K] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                  if (PREFAC) {
                    NOUT.print9997('ZPOSVXX', FACT, UPLO, N, EQUED.value, IMAT,
                        K, RESULT[K]);
                  } else {
                    NOUT.print9998(
                        'ZPOSVXX', FACT, UPLO, N, IMAT, K, RESULT[K]);
                  }
                  NFAIL++;
                }
              }
              NRUN += 7 - K1;
            }
          }
        }
      }
    }
  }

  // Print a summary of the results.

  alasvm(PATH, NOUT, NFAIL, NRUN, NERRS.value);

  // Test Error Bounds for ZGESVXX

  zebchvxx(THRESH, PATH, NOUT);
}

extension on Nout {
  void print9998(String s, String fact, String uplo, int n, int type, int test,
      double ratio) {
    println(
        ' $s, FACT=\'${fact.a1}\', UPLO=\'${uplo.a1}\', N=${n.i5}, type ${type.i1}, test(${test.i1})=${ratio.g12_5}');
  }

  void print9997(String s, String fact, String uplo, int n, String equed,
      int type, int test, double ratio) {
    println(
        ' $s, FACT=\'${fact.a1}\', UPLO=\'${uplo.a1}\', N=${n.i5}, EQUED=\'${equed.a1}\', type ${type.i1}, test(${test.i1}) =${ratio.g12_5}');
  }
}
