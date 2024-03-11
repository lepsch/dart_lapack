import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgeequ.dart';
import 'package:lapack/src/dgesv.dart';
import 'package:lapack/src/dgesvx.dart';
import 'package:lapack/src/dgesvxx.dart';
import 'package:lapack/src/dgetrf.dart';
import 'package:lapack/src/dgetri.dart';
import 'package:lapack/src/dla_gerpvgrw.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlantr.dart';
import 'package:lapack/src/dlaqge.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../matgen/dlatms.dart';
import 'aladhd.dart';
import 'alaerh.dart';
import 'alasvm.dart';
import 'common.dart';
import 'debchvxx.dart';
import 'derrvxx.dart';
import 'dget01.dart';
import 'dget02.dart';
import 'dget04.dart';
import 'dget06.dart';
import 'dget07.dart';
import 'dlarhs.dart';
import 'dlatb4.dart';
import 'xlaenv.dart';

void ddrvge(
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
  const NTYPES = 11;
  const NTESTS = 7;
  const NTRAN = 3;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS),
      BERR = Array<double>(NRHS),
      ERRBNDS_N = Matrix<double>(NRHS, 3),
      ERRBNDS_C = Matrix<double>(NRHS, 3);

  const ISEEDY = [1988, 1989, 1990, 1991];
  const TRANSS = ['N', 'T', 'C'];
  const FACTS = ['F', 'N', 'E'];
  const EQUEDS = ['N', 'R', 'C', 'B'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}GE';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

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
    final N = NVAL[IN];
    final LDA = max(N, 1);
    var XTYPE = 'N';
    final NIMAT = N <= 0 ? 1 : NTYPES;

    for (var IMAT = 1; IMAT <= NIMAT; IMAT++) {
      // Do the tests only if DOTYPE( IMAT ) is true.

      if (!DOTYPE[IMAT]) continue;

      // Skip types 5, 6, or 7 if the matrix size is too small.

      final ZEROT = IMAT >= 5 && IMAT <= 7;
      if (ZEROT && N < IMAT - 4) continue;

      // Set up parameters with DLATB4 and generate a test matrix
      // with DLATMS.

      final (:TYPE, :KL, :KU, :ANORM, :MODE, COND: CNDNUM, :DIST) =
          dlatb4(PATH, IMAT, N, N);
      var RCONDC = ONE / CNDNUM;

      srnamc.SRNAMT = 'DLATMS';
      dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
          'No packing', A.asMatrix(), LDA, WORK, INFO);

      // Check error code from DLATMS.

      if (INFO.value != 0) {
        alaerh(PATH, 'DLATMS', INFO.value, 0, ' ', N, N, -1, -1, -1, IMAT,
            NFAIL, NERRS, NOUT);
        continue;
      }

      // For types 5-7, zero one or more columns of the matrix to
      // test that INFO is returned correctly.

      final int IZERO;
      if (ZEROT) {
        if (IMAT == 5) {
          IZERO = 1;
        } else if (IMAT == 6) {
          IZERO = N;
        } else {
          IZERO = N ~/ 2 + 1;
        }
        var IOFF = (IZERO - 1) * LDA;
        if (IMAT < 7) {
          for (var I = 1; I <= N; I++) {
            A[IOFF + I] = ZERO;
          }
        } else {
          dlaset('Full', N, N - IZERO + 1, ZERO, ZERO, A(IOFF + 1).asMatrix(),
              LDA);
        }
      } else {
        IZERO = 0;
      }

      // Save a copy of the matrix A in ASAV.

      dlacpy('Full', N, N, A.asMatrix(), LDA, ASAV.asMatrix(), LDA);

      for (var IEQUED = 1; IEQUED <= 4; IEQUED++) {
        final EQUED = Box(EQUEDS[IEQUED - 1]);
        final NFACT = IEQUED == 1 ? 3 : 1;

        var RCONDO = ZERO, RCONDI = ZERO, ROLDO = ZERO, ROLDI = ZERO;
        final ROWCND = Box(ZERO), COLCND = Box(ZERO), AMAX = Box(ZERO);
        for (var IFACT = 1; IFACT <= NFACT; IFACT++) {
          final FACT = FACTS[IFACT - 1];
          final PREFAC = lsame(FACT, 'F');
          final NOFACT = lsame(FACT, 'N');
          final EQUIL = lsame(FACT, 'E');

          if (ZEROT) {
            if (PREFAC) continue;
            RCONDO = ZERO;
            RCONDI = ZERO;
          } else if (!NOFACT) {
            // Compute the condition number for comparison with
            // the value returned by DGESVX (FACT = 'N' reuses
            // the condition number from the previous iteration
            // with FACT = 'F').

            dlacpy('Full', N, N, ASAV.asMatrix(), LDA, AFAC.asMatrix(), LDA);
            if (EQUIL || IEQUED > 1) {
              // Compute row and column scale factors to
              // equilibrate the matrix A.

              dgeequ(N, N, AFAC.asMatrix(), LDA, S, S(N + 1), ROWCND, COLCND,
                  AMAX, INFO);
              if (INFO.value == 0 && N > 0) {
                if (lsame(EQUED.value, 'R')) {
                  ROWCND.value = ZERO;
                  COLCND.value = ONE;
                } else if (lsame(EQUED.value, 'C')) {
                  ROWCND.value = ONE;
                  COLCND.value = ZERO;
                } else if (lsame(EQUED.value, 'B')) {
                  ROWCND.value = ZERO;
                  COLCND.value = ZERO;
                }

                // Equilibrate the matrix.

                dlaqge(N, N, AFAC.asMatrix(), LDA, S, S(N + 1), ROWCND.value,
                    COLCND.value, AMAX.value, EQUED);
              }
            }

            // Save the condition number of the non-equilibrated
            // system for use in DGET04.

            if (EQUIL) {
              ROLDO = RCONDO;
              ROLDI = RCONDI;
            }

            // Compute the 1-norm and infinity-norm of A.

            final ANORMO = dlange('1', N, N, AFAC.asMatrix(), LDA, RWORK);
            final ANORMI = dlange('I', N, N, AFAC.asMatrix(), LDA, RWORK);

            // Factor the matrix A.

            dgetrf(N, N, AFAC.asMatrix(), LDA, IWORK, INFO);

            // Form the inverse of A.

            dlacpy('Full', N, N, AFAC.asMatrix(), LDA, A.asMatrix(), LDA);
            final LWORK = NMAX * max(3, NRHS).toInt();
            dgetri(N, A.asMatrix(), LDA, IWORK, WORK, LWORK, INFO);

            // Compute the 1-norm condition number of A.

            var AINVNM = dlange('1', N, N, A.asMatrix(), LDA, RWORK);
            if (ANORMO <= ZERO || AINVNM <= ZERO) {
              RCONDO = ONE;
            } else {
              RCONDO = (ONE / ANORMO) / AINVNM;
            }

            // Compute the infinity-norm condition number of A.

            AINVNM = dlange('I', N, N, A.asMatrix(), LDA, RWORK);
            if (ANORMI <= ZERO || AINVNM <= ZERO) {
              RCONDI = ONE;
            } else {
              RCONDI = (ONE / ANORMI) / AINVNM;
            }
          }

          for (var ITRAN = 1; ITRAN <= NTRAN; ITRAN++) {
            // Do for each value of TRANS.

            final TRANS = TRANSS[ITRAN - 1];
            if (ITRAN == 1) {
              RCONDC = RCONDO;
            } else {
              RCONDC = RCONDI;
            }

            // Restore the matrix A.

            dlacpy('Full', N, N, ASAV.asMatrix(), LDA, A.asMatrix(), LDA);

            // Form an exact solution and set the right hand side.

            srnamc.SRNAMT = 'DLARHS';
            dlarhs(PATH, XTYPE, 'Full', TRANS, N, N, KL, KU, NRHS, A.asMatrix(),
                LDA, XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
            XTYPE = 'C';
            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, BSAV.asMatrix(), LDA);

            if (NOFACT && ITRAN == 1) {
              // --- Test DGESV  ---

              // Compute the LU factorization of the matrix and
              // solve the system.

              dlacpy('Full', N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);
              dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

              srnamc.SRNAMT = 'DGESV ';
              dgesv(N, NRHS, AFAC.asMatrix(), LDA, IWORK, X.asMatrix(), LDA,
                  INFO);

              // Check error code from DGESV .

              if (INFO.value != IZERO) {
                alaerh(PATH, 'DGESV ', INFO.value, IZERO, ' ', N, N, -1, -1,
                    NRHS, IMAT, NFAIL, NERRS, NOUT);
              }

              // Reconstruct matrix from factors and compute
              // residual.

              dget01(N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA, IWORK,
                  RWORK, RESULT(1));
              final int NT;
              if (IZERO == 0) {
                // Compute residual of the computed solution.

                dlacpy(
                    'Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
                dget02('No transpose', N, N, NRHS, A.asMatrix(), LDA,
                    X.asMatrix(), LDA, WORK.asMatrix(), LDA, RWORK, RESULT(2));

                // Check solution from generated exact solution.

                dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                    RESULT(3));
                NT = 3;
              } else {
                NT = 1;
              }

              // Print information about the tests that did not
              // pass the threshold.

              for (var K = 1; K <= NT; K++) {
                if (RESULT[K] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                  NOUT.println(
                      ' DGESV , N =${N.i5}, type ${IMAT.i2}, test(${K.i2}) =${RESULT[K].g12_5}');
                  NFAIL++;
                }
              }
              NRUN += NT;
            }

            {
              // --- Test DGESVX ---

              if (!PREFAC) {
                dlaset('Full', N, N, ZERO, ZERO, AFAC.asMatrix(), LDA);
              }
              dlaset('Full', N, NRHS, ZERO, ZERO, X.asMatrix(), LDA);
              if (IEQUED > 1 && N > 0) {
                // Equilibrate the matrix if FACT = 'F' and
                // EQUED.value = 'R', 'C', or 'B'.

                dlaqge(N, N, A.asMatrix(), LDA, S, S(N + 1), ROWCND.value,
                    COLCND.value, AMAX.value, EQUED);
              }

              // Solve the system and compute the condition number
              // and error bounds using DGESVX.

              final RCOND = Box(ZERO);
              srnamc.SRNAMT = 'DGESVX';
              dgesvx(
                  FACT,
                  TRANS,
                  N,
                  NRHS,
                  A.asMatrix(),
                  LDA,
                  AFAC.asMatrix(),
                  LDA,
                  IWORK,
                  EQUED,
                  S,
                  S(N + 1),
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

              // Check the error code from DGESVX.

              if (INFO.value != IZERO) {
                alaerh(PATH, 'DGESVX', INFO.value, IZERO, FACT + TRANS, N, N,
                    -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT);
              }

              // Compare WORK(1) from DGESVX with the computed
              // reciprocal pivot growth factor RPVGRW

              double RPVGRW;
              if (INFO.value != 0) {
                RPVGRW = dlantr('M', 'U', 'N', INFO.value, INFO.value,
                    AFAC.asMatrix(), LDA, WORK);
                if (RPVGRW == ZERO) {
                  RPVGRW = ONE;
                } else {
                  RPVGRW = dlange('M', N, INFO.value, A.asMatrix(), LDA, WORK) /
                      RPVGRW;
                }
              } else {
                RPVGRW =
                    dlantr('M', 'U', 'N', N, N, AFAC.asMatrix(), LDA, WORK);
                if (RPVGRW == ZERO) {
                  RPVGRW = ONE;
                } else {
                  RPVGRW = dlange('M', N, N, A.asMatrix(), LDA, WORK) / RPVGRW;
                }
              }
              RESULT[7] =
                  (RPVGRW - WORK[1]).abs() / max(WORK[1], RPVGRW) / dlamch('E');

              final int K1;
              if (!PREFAC) {
                // Reconstruct matrix from factors and compute
                // residual.

                dget01(N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA, IWORK,
                    RWORK(2 * NRHS + 1), RESULT(1));
                K1 = 1;
              } else {
                K1 = 2;
              }

              final bool TRFCON;
              if (INFO.value == 0) {
                TRFCON = false;

                // Compute residual of the computed solution.

                dlacpy('Full', N, NRHS, BSAV.asMatrix(), LDA, WORK.asMatrix(),
                    LDA);
                dget02(TRANS, N, N, NRHS, ASAV.asMatrix(), LDA, X.asMatrix(),
                    LDA, WORK.asMatrix(), LDA, RWORK(2 * NRHS + 1), RESULT(2));

                // Check solution from generated exact solution.

                if (NOFACT || (PREFAC && lsame(EQUED.value, 'N'))) {
                  dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                      RCONDC, RESULT(3));
                } else {
                  final ROLDC = ITRAN == 1 ? ROLDO : ROLDI;
                  dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                      ROLDC, RESULT(3));
                }

                // Check the error bounds from iterative
                // refinement.

                dget07(
                    TRANS,
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
                    true,
                    RWORK(NRHS + 1),
                    RESULT(4));
              } else {
                TRFCON = true;
              }

              // Compare RCOND from DGESVX with the computed value
              // in RCONDC.

              RESULT[6] = dget06(RCOND.value, RCONDC);

              // Print information about the tests that did not pass
              // the threshold.

              if (!TRFCON) {
                for (var K = K1; K <= NTESTS; K++) {
                  if (RESULT[K] >= THRESH) {
                    if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                    if (PREFAC) {
                      NOUT.print9997('DGESVX', FACT, TRANS, N, EQUED.value,
                          IMAT, K, RESULT[K]);
                    } else {
                      NOUT.print9998(
                          'DGESVX', FACT, TRANS, N, IMAT, K, RESULT[K]);
                    }
                    NFAIL++;
                  }
                }
                NRUN += 7 - K1;
              } else {
                if (RESULT[1] >= THRESH && !PREFAC) {
                  if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                  if (PREFAC) {
                    NOUT.print9997('DGESVX', FACT, TRANS, N, EQUED.value, IMAT,
                        1, RESULT[1]);
                  } else {
                    NOUT.print9998(
                        'DGESVX', FACT, TRANS, N, IMAT, 1, RESULT[1]);
                  }
                  NFAIL++;
                  NRUN++;
                }
                if (RESULT[6] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                  if (PREFAC) {
                    NOUT.print9997('DGESVX', FACT, TRANS, N, EQUED.value, IMAT,
                        6, RESULT[6]);
                  } else {
                    NOUT.print9998(
                        'DGESVX', FACT, TRANS, N, IMAT, 6, RESULT[6]);
                  }
                  NFAIL++;
                  NRUN++;
                }
                if (RESULT[7] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                  if (PREFAC) {
                    NOUT.print9997('DGESVX', FACT, TRANS, N, EQUED.value, IMAT,
                        7, RESULT[7]);
                  } else {
                    NOUT.print9998(
                        'DGESVX', FACT, TRANS, N, IMAT, 7, RESULT[7]);
                  }
                  NFAIL++;
                  NRUN++;
                }
              }
            }

            {
              // --- Test DGESVXX ---

              // Restore the matrices A and B.

              dlacpy('Full', N, N, ASAV.asMatrix(), LDA, A.asMatrix(), LDA);
              dlacpy('Full', N, NRHS, BSAV.asMatrix(), LDA, B.asMatrix(), LDA);
              if (!PREFAC) {
                dlaset('Full', N, N, ZERO, ZERO, AFAC.asMatrix(), LDA);
              }
              dlaset('Full', N, NRHS, ZERO, ZERO, X.asMatrix(), LDA);
              if (IEQUED > 1 && N > 0) {
                // Equilibrate the matrix if FACT = 'F' and
                // EQUED.value = 'R', 'C', or 'B'.

                dlaqge(N, N, A.asMatrix(), LDA, S, S(N + 1), ROWCND.value,
                    COLCND.value, AMAX.value, EQUED);
              }

              // Solve the system and compute the condition number
              // and error bounds using DGESVXX.

              final RPVGRW_SVXX = Box(ZERO), RCOND = Box(ZERO);
              srnamc.SRNAMT = 'DGESVXX';
              const N_ERR_BNDS = 3;
              dgesvxx(
                  FACT,
                  TRANS,
                  N,
                  NRHS,
                  A.asMatrix(),
                  LDA,
                  AFAC.asMatrix(),
                  LDA,
                  IWORK,
                  EQUED,
                  S,
                  S(N + 1),
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

              // Check the error code from DGESVXX.

              if (INFO.value == N + 1) continue;
              if (INFO.value != IZERO) {
                alaerh(PATH, 'DGESVXX', INFO.value, IZERO, FACT + TRANS, N, N,
                    -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT);
                continue;
              }

              // Compare rpvgrw_svxx from DGESVXX with the computed
              // reciprocal pivot growth factor RPVGRW

              final double RPVGRW;
              if (INFO.value > 0 && INFO.value < N + 1) {
                RPVGRW = dla_gerpvgrw(
                    N, INFO.value, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);
              } else {
                RPVGRW =
                    dla_gerpvgrw(N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);
              }
              RESULT[7] = (RPVGRW - RPVGRW_SVXX.value).abs() /
                  max(RPVGRW_SVXX.value, RPVGRW) /
                  dlamch('E');

              final int K1;
              if (!PREFAC) {
                // Reconstruct matrix from factors and compute
                // residual.

                dget01(N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA, IWORK,
                    RWORK(2 * NRHS + 1), RESULT(1));
                K1 = 1;
              } else {
                K1 = 2;
              }

              final bool TRFCON;
              if (INFO.value == 0) {
                TRFCON = false;

                // Compute residual of the computed solution.

                dlacpy('Full', N, NRHS, BSAV.asMatrix(), LDA, WORK.asMatrix(),
                    LDA);
                dget02(TRANS, N, N, NRHS, ASAV.asMatrix(), LDA, X.asMatrix(),
                    LDA, WORK.asMatrix(), LDA, RWORK(2 * NRHS + 1), RESULT(2));

                // Check solution from generated exact solution.

                if (NOFACT || (PREFAC && lsame(EQUED.value, 'N'))) {
                  dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                      RCONDC, RESULT(3));
                } else {
                  final ROLDC = ITRAN == 1 ? ROLDO : ROLDI;
                  dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                      ROLDC, RESULT(3));
                }
              } else {
                TRFCON = true;
              }

              // Compare RCOND from DGESVXX with the computed value
              // in RCONDC.

              RESULT[6] = dget06(RCOND.value, RCONDC);

              // Print information about the tests that did not pass
              // the threshold.

              if (!TRFCON) {
                for (var K = K1; K <= NTESTS; K++) {
                  if (RESULT[K] >= THRESH) {
                    if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                    if (PREFAC) {
                      NOUT.print9997('DGESVXX', FACT, TRANS, N, EQUED.value,
                          IMAT, K, RESULT[K]);
                    } else {
                      NOUT.print9998(
                          'DGESVXX', FACT, TRANS, N, IMAT, K, RESULT[K]);
                    }
                    NFAIL++;
                  }
                }
                NRUN += 7 - K1;
              } else {
                if (RESULT[1] >= THRESH && !PREFAC) {
                  if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                  if (PREFAC) {
                    NOUT.print9997('DGESVXX', FACT, TRANS, N, EQUED.value, IMAT,
                        1, RESULT[1]);
                  } else {
                    NOUT.print9998(
                        'DGESVXX', FACT, TRANS, N, IMAT, 1, RESULT[1]);
                  }
                  NFAIL++;
                  NRUN++;
                }
                if (RESULT[6] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                  if (PREFAC) {
                    NOUT.print9997('DGESVXX', FACT, TRANS, N, EQUED.value, IMAT,
                        6, RESULT[6]);
                  } else {
                    NOUT.print9998(
                        'DGESVXX', FACT, TRANS, N, IMAT, 6, RESULT[6]);
                  }
                  NFAIL++;
                  NRUN++;
                }
                if (RESULT[7] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                  if (PREFAC) {
                    NOUT.print9997('DGESVXX', FACT, TRANS, N, EQUED.value, IMAT,
                        7, RESULT[7]);
                  } else {
                    NOUT.print9998(
                        'DGESVXX', FACT, TRANS, N, IMAT, 7, RESULT[7]);
                  }
                  NFAIL++;
                  NRUN++;
                }
              }
            }
          }
        }
      }
    }
  }

  // Print a summary of the results.

  alasvm(PATH, NOUT, NFAIL, NRUN, NERRS.value);

  // Test Error Bounds from DGESVXX

  debchvxx(THRESH, PATH);
}

extension on Nout {
  void print9998(String s, String fact, String trans, int n, int type, int test,
      double ratio) {
    println(
        ' $s, FACT=\'${fact.a1}\', TRANS=\'${trans.a1}\', N=${n.i5}, type ${type.i2}, test(${test.i1})=${ratio.g12_5}');
  }

  void print9997(String s, String fact, String trans, int n, String equed,
      int type, int test, double ratio) {
    println(
        ' $s, FACT=\'${fact.a1}\', TRANS=\'${trans.a1}\', N=${n.i5}, EQUED.value=\'${equed.a1}\', type ${type.i2}, test(${test.i1})=${ratio.g12_5}');
  }
}
