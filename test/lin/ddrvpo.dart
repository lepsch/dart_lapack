import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dlaqsy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dpoequ.dart';
import 'package:lapack/src/dposv.dart';
import 'package:lapack/src/dposvx.dart';
import 'package:lapack/src/dpotrf.dart';
import 'package:lapack/src/dpotri.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/lsame.dart';
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
import 'dpot01.dart';
import 'dpot02.dart';
import 'dpot05.dart';
import 'xlaenv.dart';

void ddrvpo(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NRHS,
  final double THRESH,
  final bool TSTERR,
  final int NMAX,
  final Array<double> A_,
  final Array<double> AFAC_,
  final Array<double> ASAV_,
  final Array<double> B_,
  final Array<double> BSAV_,
  final Array<double> X_,
  final Array<double> XACT_,
  final Array<double> S_,
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
  final ASAV = ASAV_.having();
  final B = B_.having();
  final BSAV = BSAV_.having();
  final X = X_.having();
  final XACT = XACT_.having();
  final S = S_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  const NTYPES = 9;
  const NTESTS = 6;
  int NFAIL, NRUN;
  double RCONDC = 0, ROLDC = 0;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'];
  const FACTS = ['F', 'N', 'E'];
  const EQUEDS = ['N', 'Y'];
  final INFO = Box(0), NERRS = Box(0);
  final SCOND = Box(0.0), AMAX = Box(0.0);

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}PO';
  NRUN = 0;
  NFAIL = 0;
  NERRS.value = 0;
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

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

      // Skip types 3, 4, or 5 if the matrix size is too small.

      final ZEROT = IMAT >= 3 && IMAT <= 5;
      if (ZEROT && N < IMAT - 2) continue;

      // Do first for UPLO = 'U', then for UPLO = 'L'

      for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
        final UPLO = UPLOS[IUPLO - 1];

        // Set up parameters with DLATB4 and generate a test matrix
        // with DLATMS.

        final (:TYPE, :KL, :KU, :ANORM, :MODE, :COND, :DIST) =
            dlatb4(PATH, IMAT, N, N);

        srnamc.SRNAMT = 'DLATMS';
        dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, COND, ANORM, KL, KU, UPLO,
            A.asMatrix(), LDA, WORK, INFO);

        // Check error code from DLATMS.

        if (INFO.value != 0) {
          alaerh(PATH, 'DLATMS', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
              NFAIL, NERRS, NOUT);
          continue;
        }

        // For types 3-5, zero one row and column of the matrix to
        // test that INFO.value is returned correctly.

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
              A[IOFF + I] = ZERO;
            }
            IOFF += IZERO;
            for (var I = IZERO; I <= N; I++) {
              A[IOFF] = ZERO;
              IOFF += LDA;
            }
          } else {
            IOFF = IZERO;
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
          IZERO = 0;
        }

        // Save a copy of the matrix A in ASAV.

        dlacpy(UPLO, N, N, A.asMatrix(), LDA, ASAV.asMatrix(), LDA);

        for (var IEQUED = 1; IEQUED <= 2; IEQUED++) {
          final EQUED = Box(EQUEDS[IEQUED - 1]);
          final NFACT = IEQUED == 1 ? 3 : 1;

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
              // the value returned by DPOSVX (FACT = 'N' reuses
              // the condition number from the previous iteration
              // with FACT = 'F').

              dlacpy(UPLO, N, N, ASAV.asMatrix(), LDA, AFAC.asMatrix(), LDA);
              if (EQUIL || IEQUED > 1) {
                // Compute row and column scale factors to
                // equilibrate the matrix A.

                dpoequ(N, AFAC.asMatrix(), LDA, S, SCOND, AMAX, INFO);
                if (INFO.value == 0 && N > 0) {
                  if (IEQUED > 1) SCOND.value = ZERO;

                  // Equilibrate the matrix.

                  dlaqsy(UPLO, N, AFAC.asMatrix(), LDA, S, SCOND.value,
                      AMAX.value, EQUED);
                }
              }

              // Save the condition number of the
              // non-equilibrated system for use in DGET04.

              if (EQUIL) ROLDC = RCONDC;

              // Compute the 1-norm of A.

              final ANORM = dlansy('1', UPLO, N, AFAC.asMatrix(), LDA, RWORK);

              // Factor the matrix A.

              dpotrf(UPLO, N, AFAC.asMatrix(), LDA, INFO);

              // Form the inverse of A.

              dlacpy(UPLO, N, N, AFAC.asMatrix(), LDA, A.asMatrix(), LDA);
              dpotri(UPLO, N, A.asMatrix(), LDA, INFO);

              // Compute the 1-norm condition number of A.

              final AINVNM = dlansy('1', UPLO, N, A.asMatrix(), LDA, RWORK);
              if (ANORM <= ZERO || AINVNM <= ZERO) {
                RCONDC = ONE;
              } else {
                RCONDC = (ONE / ANORM) / AINVNM;
              }
            }

            // Restore the matrix A.

            dlacpy(UPLO, N, N, ASAV.asMatrix(), LDA, A.asMatrix(), LDA);

            // Form an exact solution and set the right hand side.

            srnamc.SRNAMT = 'DLARHS';
            dlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A.asMatrix(),
                LDA, XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
            XTYPE = 'C';
            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, BSAV.asMatrix(), LDA);

            if (NOFACT) {
              // --- Test DPOSV  ---

              // Compute the L*L' or U'*U factorization of the
              // matrix and solve the system.

              dlacpy(UPLO, N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);
              dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

              srnamc.SRNAMT = 'DPOSV ';
              dposv(
                  UPLO, N, NRHS, AFAC.asMatrix(), LDA, X.asMatrix(), LDA, INFO);

              // Check error code from DPOSV .

              if (INFO.value != IZERO) {
                alaerh(PATH, 'DPOSV ', INFO.value, IZERO, UPLO, N, N, -1, -1,
                    NRHS, IMAT, NFAIL, NERRS, NOUT);
                //
              } else if (INFO.value != 0) {
                //
              } else {
                // Reconstruct matrix from factors and compute
                // residual.

                dpot01(UPLO, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA, RWORK,
                    RESULT(1));

                // Compute residual of the computed solution.

                dlacpy(
                    'Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
                dpot02(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                    WORK.asMatrix(), LDA, RWORK, RESULT(2));

                // Check solution from generated exact solution.

                dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                    RESULT(3));
                final NT = 3;

                // Print information about the tests that did not
                // pass the threshold.

                for (var K = 1; K <= NT; K++) {
                  if (RESULT[K] >= THRESH) {
                    if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                    NOUT.println(
                        ' DPOSV , UPLO=\'${UPLO.a1}\', N =${N.i5}, type ${IMAT.i1}, test(${K.i1})=${RESULT[K].g12_5}');
                    NFAIL++;
                  }
                }
                NRUN += NT;
              }
            }

            // --- Test DPOSVX ---

            if (!PREFAC) dlaset(UPLO, N, N, ZERO, ZERO, AFAC.asMatrix(), LDA);
            dlaset('Full', N, NRHS, ZERO, ZERO, X.asMatrix(), LDA);
            if (IEQUED > 1 && N > 0) {
              // Equilibrate the matrix if FACT='F' and
              // EQUED.value='Y'.

              dlaqsy(UPLO, N, A.asMatrix(), LDA, S, SCOND.value, AMAX.value,
                  EQUED);
            }

            // Solve the system and compute the condition number
            // and error bounds using DPOSVX.

            final RCOND = Box(0.0);
            srnamc.SRNAMT = 'DPOSVX';
            dposvx(
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
                IWORK,
                INFO);

            // Check the error code from DPOSVX.

            if (INFO.value != IZERO) {
              alaerh(PATH, 'DPOSVX', INFO.value, IZERO, FACT + UPLO, N, N, -1,
                  -1, NRHS, IMAT, NFAIL, NERRS, NOUT);
              continue;
            }

            final int K1;
            if (INFO.value == 0) {
              if (!PREFAC) {
                // Reconstruct matrix from factors and compute
                // residual.

                dpot01(UPLO, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA,
                    RWORK(2 * NRHS + 1), RESULT(1));
                K1 = 1;
              } else {
                K1 = 2;
              }

              // Compute residual of the computed solution.

              dlacpy(
                  'Full', N, NRHS, BSAV.asMatrix(), LDA, WORK.asMatrix(), LDA);
              dpot02(UPLO, N, NRHS, ASAV.asMatrix(), LDA, X.asMatrix(), LDA,
                  WORK.asMatrix(), LDA, RWORK(2 * NRHS + 1), RESULT(2));

              // Check solution from generated exact solution.

              if (NOFACT || (PREFAC && lsame(EQUED.value, 'N'))) {
                dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                    RESULT(3));
              } else {
                dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, ROLDC,
                    RESULT(3));
              }

              // Check the error bounds from iterative
              // refinement.

              dpot05(
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

            // Compare RCOND from DPOSVX with the computed value
            // in RCONDC.

            RESULT[6] = dget06(RCOND.value, RCONDC);

            // Print information about the tests that did not pass
            // the threshold.

            for (var K = K1; K <= 6; K++) {
              if (RESULT[K] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                if (PREFAC) {
                  NOUT.println(
                      ' DPOSVX, FACT=\'${FACT.a1}\', UPLO=\'${UPLO.a1}\', N=${N.i5}, EQUED=\'${EQUED.value.a1}\', type ${IMAT.i1}, test(${K.i1}) =${RESULT[K].g12_5}');
                } else {
                  NOUT.println(
                      ' DPOSVX, FACT=\'${FACT.a1}\', UPLO=\'${UPLO.a1}\', N=${N.i5}, type ${IMAT.i1}, test(${K.i1})=${RESULT[K].g12_5}');
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

  // Print a summary of the results.

  alasvm(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
