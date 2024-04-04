import 'dart:math';

import 'package:lapack/lapack.dart';

import '../matgen/zlatms.dart';
import 'aladhd.dart';
import 'alaerh.dart';
import 'alasvm.dart';
import 'common.dart';
import 'dget06.dart';
import 'zerrvx.dart';
import 'zget04.dart';
import 'zlaipd.dart';
import 'zlarhs.dart';
import 'zlatb4.dart';
import 'zppt01.dart';
import 'zppt02.dart';
import 'zppt05.dart';

void zdrvpp(
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
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'],
      FACTS = ['F', 'N', 'E'],
      PACKS = ['C', 'R'],
      EQUEDS = ['N', 'Y'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Zomplex precision'[0]}PP';
  var NRUN = 0;
  var NFAIL = 0;
  var NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) zerrvx(PATH, NOUT);
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

      // Skip types 3, 4, or 5 if the matrix size is too small.

      final ZEROT = IMAT >= 3 && IMAT <= 5;
      if (ZEROT && N < IMAT - 2) continue;

      // Do first for UPLO = 'U', then for UPLO = 'L'

      for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
        final UPLO = UPLOS[IUPLO - 1];
        final PACKIT = PACKS[IUPLO - 1];

        // Set up parameters with ZLATB4 and generate a test matrix
        // with ZLATMS.

        final (:TYPE, :KL, :KU, :ANORM, :MODE, :CNDNUM, :DIST) =
            zlatb4(PATH, IMAT, N, N);
        var RCONDC = ONE / CNDNUM;

        srnamc.SRNAMT = 'ZLATMS';
        zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
            PACKIT, A.asMatrix(), LDA, WORK, INFO);

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

          // Set row and column IZERO of A to 0.

          if (IUPLO == 1) {
            var IOFF = (IZERO - 1) * IZERO ~/ 2;
            for (var I = 1; I <= IZERO - 1; I++) {
              A[IOFF + I] = Complex.zero;
            }
            IOFF += IZERO;
            for (var I = IZERO; I <= N; I++) {
              A[IOFF] = Complex.zero;
              IOFF += I;
            }
          } else {
            var IOFF = IZERO;
            for (var I = 1; I <= IZERO - 1; I++) {
              A[IOFF] = Complex.zero;
              IOFF += N - I;
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

        if (IUPLO == 1) {
          zlaipd(N, A, 2, 1);
        } else {
          zlaipd(N, A, N, -1);
        }

        // Save a copy of the matrix A in ASAV.

        zcopy(NPP, A, 1, ASAV, 1);

        for (var IEQUED = 1; IEQUED <= 2; IEQUED++) {
          final EQUED = Box(EQUEDS[IEQUED - 1]);
          final NFACT = IEQUED == 1 ? 3 : 1;

          double ROLDC = ZERO;
          for (var IFACT = 1; IFACT <= NFACT; IFACT++) {
            final FACT = FACTS[IFACT - 1];
            final PREFAC = lsame(FACT, 'F');
            final NOFACT = lsame(FACT, 'N');
            final EQUIL = lsame(FACT, 'E');

            final SCOND = Box(ZERO), AMAX = Box(ZERO);
            if (ZEROT) {
              if (PREFAC) continue;
              RCONDC = ZERO;
            } else if (!lsame(FACT, 'N')) {
              // Compute the condition number for comparison with
              // the value returned by ZPPSVX (FACT = 'N' reuses
              // the condition number from the previous iteration
              //    with FACT = 'F').

              zcopy(NPP, ASAV, 1, AFAC, 1);
              if (EQUIL || IEQUED > 1) {
                // Compute row and column scale factors to
                // equilibrate the matrix A.

                zppequ(UPLO, N, AFAC, S, SCOND, AMAX, INFO);
                if (INFO.value == 0 && N > 0) {
                  if (IEQUED > 1) SCOND.value = ZERO;

                  // Equilibrate the matrix.

                  zlaqhp(UPLO, N, AFAC, S, SCOND.value, AMAX.value, EQUED);
                }
              }

              // Save the condition number of the
              // non-equilibrated system for use in ZGET04.

              if (EQUIL) ROLDC = RCONDC;

              // Compute the 1-norm of A.

              final ANORM = zlanhp('1', UPLO, N, AFAC, RWORK);

              // Factor the matrix A.

              zpptrf(UPLO, N, AFAC, INFO);

              // Form the inverse of A.

              zcopy(NPP, AFAC, 1, A, 1);
              zpptri(UPLO, N, A, INFO);

              // Compute the 1-norm condition number of A.

              final AINVNM = zlanhp('1', UPLO, N, A, RWORK);
              if (ANORM <= ZERO || AINVNM <= ZERO) {
                RCONDC = ONE;
              } else {
                RCONDC = (ONE / ANORM) / AINVNM;
              }
            }

            // Restore the matrix A.

            zcopy(NPP, ASAV, 1, A, 1);

            // Form an exact solution and set the right hand side.

            srnamc.SRNAMT = 'ZLARHS';
            zlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A.asMatrix(),
                LDA, XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
            XTYPE = 'C';
            zlacpy('Full', N, NRHS, B.asMatrix(), LDA, BSAV.asMatrix(), LDA);

            if (NOFACT) {
              // --- Test ZPPSV  ---

              // Compute the L*L' or U'*U factorization of the
              // matrix and solve the system.

              zcopy(NPP, A, 1, AFAC, 1);
              zlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

              srnamc.SRNAMT = 'ZPPSV';
              zppsv(UPLO, N, NRHS, AFAC, X.asMatrix(), LDA, INFO);

              // Check error code from ZPPSV .

              if (INFO.value != IZERO) {
                alaerh(PATH, 'ZPPSV ', INFO.value, IZERO, UPLO, N, N, -1, -1,
                    NRHS, IMAT, NFAIL, NERRS, NOUT);
              } else if (INFO.value != 0) {
                //
              } else {
                // Reconstruct matrix from factors and compute
                // residual.

                zppt01(UPLO, N, A, AFAC, RWORK, RESULT(1));

                // Compute residual of the computed solution.

                zlacpy(
                    'Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
                zppt02(UPLO, N, NRHS, A, X.asMatrix(), LDA, WORK.asMatrix(),
                    LDA, RWORK, RESULT(2));

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
                        ' ZPPSV , UPLO=\'${UPLO.a1}\', N =${N.i5}, type ${IMAT.i1}, test(${K.i1})=${RESULT[K].g12_5}');
                    NFAIL++;
                  }
                }
                NRUN += NT;
              }
            }

            // --- Test ZPPSVX ---

            if (!PREFAC && NPP > 0) {
              zlaset('Full', NPP, 1, Complex.zero, Complex.zero,
                  AFAC.asMatrix(), NPP);
            }
            zlaset(
                'Full', N, NRHS, Complex.zero, Complex.zero, X.asMatrix(), LDA);
            if (IEQUED > 1 && N > 0) {
              // Equilibrate the matrix if FACT='F' and
              // EQUED='Y'.

              zlaqhp(UPLO, N, A, S, SCOND.value, AMAX.value, EQUED);
            }

            // Solve the system and compute the condition number
            // and error bounds using ZPPSVX.

            final RCOND = Box(ZERO);
            srnamc.SRNAMT = 'ZPPSVX';
            zppsvx(
                FACT,
                UPLO,
                N,
                NRHS,
                A,
                AFAC,
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

            // Check the error code from ZPPSVX.

            if (INFO.value != IZERO) {
              alaerh(PATH, 'ZPPSVX', INFO.value, IZERO, FACT + UPLO, N, N, -1,
                  -1, NRHS, IMAT, NFAIL, NERRS, NOUT);
            } else {
              final int K1;
              if (INFO.value == 0) {
                if (!PREFAC) {
                  // Reconstruct matrix from factors and compute
                  // residual.

                  zppt01(UPLO, N, A, AFAC, RWORK(2 * NRHS + 1), RESULT(1));
                  K1 = 1;
                } else {
                  K1 = 2;
                }

                // Compute residual of the computed solution.

                zlacpy('Full', N, NRHS, BSAV.asMatrix(), LDA, WORK.asMatrix(),
                    LDA);
                zppt02(UPLO, N, NRHS, ASAV, X.asMatrix(), LDA, WORK.asMatrix(),
                    LDA, RWORK(2 * NRHS + 1), RESULT(2));

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

                zppt05(
                    UPLO,
                    N,
                    NRHS,
                    ASAV,
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

              // Compare RCOND from ZPPSVX with the computed value
              // in RCONDC.

              RESULT[6] = dget06(RCOND.value, RCONDC);

              // Print information about the tests that did not pass
              // the threshold.

              for (var K = K1; K <= 6; K++) {
                if (RESULT[K] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                  if (PREFAC) {
                    NOUT.println(
                        ' ZPPSVX, FACT=\'${FACT.a1}\', UPLO=\'${UPLO.a1}\', N=${N.i5}, EQUED=\'${EQUED.value.a1}\', type ${IMAT.i1}, test(${K.i1})=${RESULT[K].g12_5}');
                  } else {
                    NOUT.println(
                        ' ZPPSVX, FACT=\'${FACT.a1}\', UPLO=\'${UPLO.a1}\', N=${N.i5}, type ${IMAT.i1}, test(${K.i1})=${RESULT[K].g12_5}');
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
}
