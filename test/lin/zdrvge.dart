import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zgetrf.dart';
import 'package:lapack/src/zgeequ.dart';
import 'package:lapack/src/zgesv.dart';
import 'package:lapack/src/zgesvx.dart';
import 'package:lapack/src/zgetri.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlantr.dart';
import 'package:lapack/src/zlaqge.dart';
import 'package:lapack/src/zlaset.dart';

import '../matgen/zlatms.dart';
import 'aladhd.dart';
import 'alaerh.dart';
import 'alasvm.dart';
import 'common.dart';
import 'dget06.dart';
import 'xlaenv.dart';
import 'zerrvxx.dart';
import 'zget01.dart';
import 'zget02.dart';
import 'zget04.dart';
import 'zget07.dart';
import 'zlarhs.dart';
import 'zlatb4.dart';

void zdrvge(
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
  final RDUM = Array<double>(1), RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const TRANSS = ['N', 'T', 'C'];
  const FACTS = ['F', 'N', 'E'];
  const EQUEDS = ['N', 'R', 'C', 'B'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Zomplex precision'[0]}GE';
  var NRUN = 0;
  var NFAIL = 0;
  var NERRS = Box(0);
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

      // Skip types 5, 6, or 7 if the matrix size is too small.

      final ZEROT = IMAT >= 5 && IMAT <= 7;
      if (ZEROT && N < IMAT - 4) continue;

      // Set up parameters with ZLATB4 and generate a test matrix
      // with ZLATMS.

      final (:TYPE, :KL, :KU, :ANORM, :MODE, :CNDNUM, :DIST) =
          zlatb4(PATH, IMAT, N, N);
      var RCONDC = ONE / CNDNUM;

      srnamc.SRNAMT = 'ZLATMS';
      zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
          'No packing', A.asMatrix(), LDA, WORK, INFO);

      // Check error code from ZLATMS.

      if (INFO.value != 0) {
        alaerh(PATH, 'ZLATMS', INFO.value, 0, ' ', N, N, -1, -1, -1, IMAT,
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
            A[IOFF + I] = Complex.zero;
          }
        } else {
          zlaset('Full', N, N - IZERO + 1, Complex.zero, Complex.zero,
              A(IOFF + 1).asMatrix(), LDA);
        }
      } else {
        IZERO = 0;
      }

      // Save a copy of the matrix A in ASAV.

      zlacpy('Full', N, N, A.asMatrix(), LDA, ASAV.asMatrix(), LDA);

      for (var IEQUED = 1; IEQUED <= 4; IEQUED++) {
        final EQUED = Box(EQUEDS[IEQUED - 1]);
        final NFACT = IEQUED == 1 ? 3 : 1;

        double RCONDO = ZERO, RCONDI = ZERO;
        for (var IFACT = 1; IFACT <= NFACT; IFACT++) {
          final FACT = FACTS[IFACT - 1];
          final PREFAC = lsame(FACT, 'F');
          final NOFACT = lsame(FACT, 'N');
          final EQUIL = lsame(FACT, 'E');

          double ROLDO = ZERO, ROLDI = ZERO;
          final ROWCND = Box(ZERO), COLCND = Box(ZERO), AMAX = Box(ZERO);
          if (ZEROT) {
            if (PREFAC) continue;
            RCONDO = ZERO;
            RCONDI = ZERO;
          } else if (!NOFACT) {
            // Compute the condition number for comparison with
            // the value returned by ZGESVX (FACT = 'N' reuses
            // the condition number from the previous iteration
            // with FACT = 'F').

            zlacpy('Full', N, N, ASAV.asMatrix(), LDA, AFAC.asMatrix(), LDA);
            if (EQUIL || IEQUED > 1) {
              // Compute row and column scale factors to
              // equilibrate the matrix A.

              zgeequ(N, N, AFAC.asMatrix(), LDA, S, S(N + 1), ROWCND, COLCND,
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

                zlaqge(N, N, AFAC.asMatrix(), LDA, S, S(N + 1), ROWCND.value,
                    COLCND.value, AMAX.value, EQUED);
              }
            }

            // Save the condition number of the non-equilibrated
            // system for use in ZGET04.

            if (EQUIL) {
              ROLDO = RCONDO;
              ROLDI = RCONDI;
            }

            // Compute the 1-norm and infinity-norm of A.

            final ANORMO = zlange('1', N, N, AFAC.asMatrix(), LDA, RWORK);
            final ANORMI = zlange('I', N, N, AFAC.asMatrix(), LDA, RWORK);

            // Factor the matrix A.

            srnamc.SRNAMT = 'ZGETRF';
            zgetrf(N, N, AFAC.asMatrix(), LDA, IWORK, INFO);

            // Form the inverse of A.

            zlacpy('Full', N, N, AFAC.asMatrix(), LDA, A.asMatrix(), LDA);
            final LWORK = NMAX * max(3, NRHS).toInt();
            srnamc.SRNAMT = 'ZGETRI';
            zgetri(N, A.asMatrix(), LDA, IWORK, WORK, LWORK, INFO);

            // Compute the 1-norm condition number of A.

            var AINVNM = zlange('1', N, N, A.asMatrix(), LDA, RWORK);
            if (ANORMO <= ZERO || AINVNM <= ZERO) {
              RCONDO = ONE;
            } else {
              RCONDO = (ONE / ANORMO) / AINVNM;
            }

            // Compute the infinity-norm condition number of A.

            AINVNM = zlange('I', N, N, A.asMatrix(), LDA, RWORK);
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

            zlacpy('Full', N, N, ASAV.asMatrix(), LDA, A.asMatrix(), LDA);

            // Form an exact solution and set the right hand side.

            srnamc.SRNAMT = 'ZLARHS';
            zlarhs(PATH, XTYPE, 'Full', TRANS, N, N, KL, KU, NRHS, A.asMatrix(),
                LDA, XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
            XTYPE = 'C';
            zlacpy('Full', N, NRHS, B.asMatrix(), LDA, BSAV.asMatrix(), LDA);

            if (NOFACT && ITRAN == 1) {
              // --- Test ZGESV  ---

              // Compute the LU factorization of the matrix and
              // solve the system.

              zlacpy('Full', N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);
              zlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

              srnamc.SRNAMT = 'ZGESV ';
              zgesv(N, NRHS, AFAC.asMatrix(), LDA, IWORK, X.asMatrix(), LDA,
                  INFO);

              // Check error code from ZGESV .

              if (INFO.value != IZERO) {
                alaerh(PATH, 'ZGESV ', INFO.value, IZERO, ' ', N, N, -1, -1,
                    NRHS, IMAT, NFAIL, NERRS, NOUT);
              }

              // Reconstruct matrix from factors and compute
              // residual.

              zget01(N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA, IWORK,
                  RWORK, RESULT(1));
              final int NT;
              if (IZERO == 0) {
                // Compute residual of the computed solution.

                zlacpy(
                    'Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
                zget02('No transpose', N, N, NRHS, A.asMatrix(), LDA,
                    X.asMatrix(), LDA, WORK.asMatrix(), LDA, RWORK, RESULT(2));

                // Check solution from generated exact solution.

                zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
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
                      ' ZGESV , N =${N.i5}, type ${IMAT.i2}, test(${K.i2}) =${RESULT[K].g12_5}');
                  NFAIL++;
                }
              }
              NRUN += NT;
            }

            // --- Test ZGESVX ---

            if (!PREFAC) {
              zlaset('Full', N, N, Complex.zero, Complex.zero, AFAC.asMatrix(),
                  LDA);
            }
            zlaset(
                'Full', N, NRHS, Complex.zero, Complex.zero, X.asMatrix(), LDA);
            if (IEQUED > 1 && N > 0) {
              // Equilibrate the matrix if FACT = 'F' and
              // EQUED = 'R', 'C', or 'B'.

              zlaqge(N, N, A.asMatrix(), LDA, S, S(N + 1), ROWCND.value,
                  COLCND.value, AMAX.value, EQUED);
            }

            // Solve the system and compute the condition number
            // and error bounds using ZGESVX.

            final RCOND = Box(ZERO);
            srnamc.SRNAMT = 'ZGESVX';
            zgesvx(
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
                RWORK(2 * NRHS + 1),
                INFO);

            // Check the error code from ZGESVX.

            if (INFO.value != IZERO) {
              alaerh(PATH, 'ZGESVX', INFO.value, IZERO, FACT + TRANS, N, N, -1,
                  -1, NRHS, IMAT, NFAIL, NERRS, NOUT);
            }

            // Compare RWORK(2*NRHS+1) from ZGESVX with the
            // computed reciprocal pivot growth factor RPVGRW

            double RPVGRW;
            if (INFO.value != 0 && INFO.value <= N) {
              RPVGRW = zlantr('M', 'U', 'N', INFO.value, INFO.value,
                  AFAC.asMatrix(), LDA, RDUM);
              if (RPVGRW == ZERO) {
                RPVGRW = ONE;
              } else {
                RPVGRW = zlange('M', N, INFO.value, A.asMatrix(), LDA, RDUM) /
                    RPVGRW;
              }
            } else {
              RPVGRW = zlantr('M', 'U', 'N', N, N, AFAC.asMatrix(), LDA, RDUM);
              if (RPVGRW == ZERO) {
                RPVGRW = ONE;
              } else {
                RPVGRW = zlange('M', N, N, A.asMatrix(), LDA, RDUM) / RPVGRW;
              }
            }
            RESULT[7] = (RPVGRW - RWORK[2 * NRHS + 1]).abs() /
                max(RWORK[2 * NRHS + 1], RPVGRW) /
                dlamch('E');

            final int K1;
            if (!PREFAC) {
              // Reconstruct matrix from factors and compute
              // residual.

              zget01(N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA, IWORK,
                  RWORK(2 * NRHS + 1), RESULT(1));
              K1 = 1;
            } else {
              K1 = 2;
            }

            final bool TRFCON;
            if (INFO.value == 0) {
              TRFCON = false;

              // Compute residual of the computed solution.

              zlacpy(
                  'Full', N, NRHS, BSAV.asMatrix(), LDA, WORK.asMatrix(), LDA);
              zget02(TRANS, N, N, NRHS, ASAV.asMatrix(), LDA, X.asMatrix(), LDA,
                  WORK.asMatrix(), LDA, RWORK(2 * NRHS + 1), RESULT(2));

              // Check solution from generated exact solution.

              if (NOFACT || (PREFAC && lsame(EQUED.value, 'N'))) {
                zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                    RESULT(3));
              } else {
                final ROLDC = ITRAN == 1 ? ROLDO : ROLDI;
                zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, ROLDC,
                    RESULT(3));
              }

              // Check the error bounds from iterative
              // refinement.

              zget07(
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

            // Compare RCOND from ZGESVX with the computed value
            // in RCONDC.

            RESULT[6] = dget06(RCOND.value, RCONDC);

            // Print information about the tests that did not pass
            // the threshold.

            if (!TRFCON) {
              for (var K = K1; K <= NTESTS; K++) {
                if (RESULT[K] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                  if (PREFAC) {
                    NOUT.print9997('ZGESVX', FACT, TRANS, N, EQUED.value, IMAT,
                        K, RESULT[K]);
                  } else {
                    NOUT.print9998(
                        'ZGESVX', FACT, TRANS, N, IMAT, K, RESULT[K]);
                  }
                  NFAIL++;
                }
              }
              NRUN += NTESTS - K1 + 1;
            } else {
              if (RESULT[1] >= THRESH && !PREFAC) {
                if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                if (PREFAC) {
                  NOUT.print9997('ZGESVX', FACT, TRANS, N, EQUED.value, IMAT, 1,
                      RESULT[1]);
                } else {
                  NOUT.print9998('ZGESVX', FACT, TRANS, N, IMAT, 1, RESULT[1]);
                }
                NFAIL++;
                NRUN++;
              }
              if (RESULT[6] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                if (PREFAC) {
                  NOUT.print9997('ZGESVX', FACT, TRANS, N, EQUED.value, IMAT, 6,
                      RESULT[6]);
                } else {
                  NOUT.print9998('ZGESVX', FACT, TRANS, N, IMAT, 6, RESULT[6]);
                }
                NFAIL++;
                NRUN++;
              }
              if (RESULT[7] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                if (PREFAC) {
                  NOUT.print9997('ZGESVX', FACT, TRANS, N, EQUED.value, IMAT, 7,
                      RESULT[7]);
                } else {
                  NOUT.print9998('ZGESVX', FACT, TRANS, N, IMAT, 7, RESULT[7]);
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

  // Print a summary of the results.

  alasvm(PATH, NOUT, NFAIL, NRUN, NERRS.value);
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
        ' $s, FACT=\'${fact.a1}\', TRANS=\'${trans.a1}\', N=${n.i5}, EQUED=\'${equed.a1}\', type ${type.i2}, test(${test.i1})=${ratio.g12_5}');
  }
}
