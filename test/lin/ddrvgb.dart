import 'dart:math';

import 'package:lapack/lapack.dart';
import 'package:test/test.dart';

import '../matgen/dlatms.dart';
import '../test_driver.dart';
import 'aladhd.dart';
import 'alaerh.dart';
import 'alasvm.dart';
import 'common.dart';
import 'derrvx.dart';
import 'dgbt01.dart';
import 'dgbt02.dart';
import 'dgbt05.dart';
import 'dget04.dart';
import 'dget06.dart';
import 'dlarhs.dart';
import 'dlatb4.dart';
import 'xlaenv.dart';

void ddrvgb(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NRHS,
  final double THRESH,
  final bool TSTERR,
  final Array<double> A_,
  final int LA,
  final Array<double> AFB_,
  final int LAFB,
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
  final TestDriver test,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DOTYPE = DOTYPE_.having();
  final NVAL = NVAL_.having();
  final A = A_.having();
  final AFB = AFB_.having();
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
  const NTYPES = 8, NTESTS = 7, NTRAN = 3;
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const TRANSS = ['N', 'T', 'C'];
  const FACTS = ['F', 'N', 'E'];
  const EQUEDS = ['N', 'R', 'C', 'B'];

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}GB';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  final ISEED = Array.fromList(ISEEDY);

  // Test the error exits
  test.group('error exits', () {
    if (TSTERR) derrvx(PATH, NOUT, test);
    test.tearDown(() {
      infoc.INFOT = 0;
    });
  });

  test.setUp(() {
    // Set the block size and minimum block size for testing.
    const NB = 1;
    const NBMIN = 2;
    xlaenv(1, NB);
    xlaenv(2, NBMIN);
  });

  // Do for each value of N in NVAL
  for (final IN in 1.through(NN)) {
    final N = NVAL[IN];
    final LDB = max(N, 1);

    // Set limits on the number of loop iterations.
    final NKL = N == 0 ? 1 : max(1, min(N, 4));
    final NKU = NKL;
    final NIMAT = N <= 0 ? 1 : NTYPES;

    for (final IKL in 1.through(NKL)) {
      // Do for KL = 0, N-1, (3N-1)/4, and (N+1)/4. This order makes
      // it easier to skip redundant values for small values of N.

      final KL = switch (IKL) {
        1 => 0,
        2 => max(N - 1, 0),
        3 => (3 * N - 1) ~/ 4,
        4 => (N + 1) ~/ 4,
        _ => throw UnimplementedError(),
      };

      for (final IKU in 1.through(NKU)) {
        // Do for KU = 0, N-1, (3N-1)/4, and (N+1)/4. This order
        // makes it easier to skip redundant values for small
        // values of N.

        final KU = switch (IKU) {
          1 => 0,
          2 => max(N - 1, 0),
          3 => (3 * N - 1) ~/ 4,
          4 => (N + 1) ~/ 4,
          _ => throw UnimplementedError(),
        };

        // Check that A and AFB are big enough to generate this
        // matrix.

        final LDA = KL + KU + 1;
        final LDAFB = 2 * KL + KU + 1;
        if (LDA * N > LA || LDAFB * N > LAFB) {
          if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
          if (LDA * N > LA) {
            NOUT.println(
                ' *** In DDRVGB, LA=${LA.i5} is too small for N=${N.i5}, KU=${KU.i5}, KL=${KL.i5}\n ==> Increase LA to at least ${(N * (KL + KU + 1)).i5}');
            NERRS.value++;
          }
          if (LDAFB * N > LAFB) {
            NOUT.println(
                ' *** In DDRVGB, LAFB=${LAFB.i5} is too small for N=${N.i5}, KU=${KU.i5}, KL=${KL.i5}\n ==> Increase LAFB to at least ${(N * (2 * KL + KU + 1)).i5}');
            NERRS.value++;
          }
          continue;
        }

        for (final IMAT in 1.through(NIMAT)) {
          // Do the tests only if DOTYPE( IMAT ) is true.
          final skip = !DOTYPE[IMAT];

          // Skip types 2, 3, or 4 if the matrix is too small.
          final ZEROT = IMAT >= 2 && IMAT <= 4;
          if (ZEROT && N < IMAT - 1) continue;

          test('DDRVGB (IN=$IN IKL=$IKL IKU=$IKU IMAT=$IMAT)', () {
            final INFO = Box(0);
            String? XTYPE;

            // Set up parameters with DLATB4 and generate a
            // test matrix with DLATMS.

            final (:TYPE, :KL, :KU, :ANORM, :MODE, COND: CNDNUM, :DIST) =
                dlatb4(PATH, IMAT, N, N);
            var RCONDC = ONE / CNDNUM;

            srnamc.SRNAMT = 'DLATMS';
            dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
                'Z', A.asMatrix(), LDA, WORK, INFO);

            // Check the error code from DLATMS.
            test.expect(INFO.value, 0);
            if (INFO.value != 0) {
              alaerh(PATH, 'DLATMS', INFO.value, 0, ' ', N, N, KL, KU, -1, IMAT,
                  NFAIL, NERRS, NOUT);
              return;
            }

            // For types 2, 3, and 4, zero one or more columns of
            // the matrix to test that INFO is returned correctly.

            final int IZERO;
            if (ZEROT) {
              if (IMAT == 2) {
                IZERO = 1;
              } else if (IMAT == 3) {
                IZERO = N;
              } else {
                IZERO = N ~/ 2 + 1;
              }
              var IOFF = (IZERO - 1) * LDA;
              if (IMAT < 4) {
                final I1 = max(1, KU + 2 - IZERO);
                final I2 = min(KL + KU + 1, KU + 1 + (N - IZERO));
                for (var I = I1; I <= I2; I++) {
                  A[IOFF + I] = ZERO;
                }
              } else {
                for (var J = IZERO; J <= N; J++) {
                  for (var I = max(1, KU + 2 - J);
                      I <= min(KL + KU + 1, KU + 1 + (N - J));
                      I++) {
                    A[IOFF + I] = ZERO;
                  }
                  IOFF += LDA;
                }
              }
            } else {
              IZERO = 0;
            }

            // Save a copy of the matrix A in ASAV.
            dlacpy('Full', KL + KU + 1, N, A.asMatrix(), LDA, ASAV.asMatrix(),
                LDA);

            for (var IEQUED = 1; IEQUED <= 4; IEQUED++) {
              final EQUED = Box(EQUEDS[IEQUED - 1]);
              final NFACT = IEQUED == 1 ? 3 : 1;

              var RCONDO = ZERO, RCONDI = ZERO, ROLDO = ZERO, ROLDI = ZERO;
              for (var IFACT = 1; IFACT <= NFACT; IFACT++) {
                final FACT = FACTS[IFACT - 1];
                final PREFAC = lsame(FACT, 'F');
                final NOFACT = lsame(FACT, 'N');
                final EQUIL = lsame(FACT, 'E');

                final ROWCND = Box(0.0), COLCND = Box(0.0), AMAX = Box(0.0);
                if (ZEROT) {
                  if (PREFAC) continue;
                  RCONDO = ZERO;
                  RCONDI = ZERO;
                } else if (!NOFACT) {
                  // Compute the condition number for comparison
                  // with the value returned by DGESVX (FACT =
                  // 'N' reuses the condition number from the
                  // previous iteration with FACT = 'F').

                  dlacpy('Full', KL + KU + 1, N, ASAV.asMatrix(), LDA,
                      AFB(KL + 1).asMatrix(), LDAFB);
                  if (EQUIL || IEQUED > 1) {
                    // Compute row and column scale factors to
                    // equilibrate the matrix A.

                    dgbequ(N, N, KL, KU, AFB(KL + 1).asMatrix(), LDAFB, S,
                        S(N + 1), ROWCND, COLCND, AMAX, INFO);
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

                      dlaqgb(
                          N,
                          N,
                          KL,
                          KU,
                          AFB(KL + 1).asMatrix(),
                          LDAFB,
                          S,
                          S(N + 1),
                          ROWCND.value,
                          COLCND.value,
                          AMAX.value,
                          EQUED);
                    }
                  }

                  // Save the condition number of the
                  // non-equilibrated system for use in DGET04.

                  if (EQUIL) {
                    ROLDO = RCONDO;
                    ROLDI = RCONDI;
                  }

                  // Compute the 1-norm and infinity-norm of A.

                  final ANORMO = dlangb(
                      '1', N, KL, KU, AFB(KL + 1).asMatrix(), LDAFB, RWORK);
                  final ANORMI = dlangb(
                      'I', N, KL, KU, AFB(KL + 1).asMatrix(), LDAFB, RWORK);

                  // Factor the matrix A.

                  dgbtrf(N, N, KL, KU, AFB.asMatrix(), LDAFB, IWORK, INFO);

                  // Form the inverse of A.

                  dlaset('Full', N, N, ZERO, ONE, WORK.asMatrix(), LDB);
                  srnamc.SRNAMT = 'DGBTRS';
                  dgbtrs('No transpose', N, KL, KU, N, AFB.asMatrix(), LDAFB,
                      IWORK, WORK.asMatrix(), LDB, INFO);

                  // Compute the 1-norm condition number of A.

                  var AINVNM = dlange('1', N, N, WORK.asMatrix(), LDB, RWORK);
                  if (ANORMO <= ZERO || AINVNM <= ZERO) {
                    RCONDO = ONE;
                  } else {
                    RCONDO = (ONE / ANORMO) / AINVNM;
                  }

                  // Compute the infinity-norm condition number
                  // of A.

                  AINVNM = dlange('I', N, N, WORK.asMatrix(), LDB, RWORK);
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

                  dlacpy('Full', KL + KU + 1, N, ASAV.asMatrix(), LDA,
                      A.asMatrix(), LDA);

                  // Form an exact solution and set the right hand
                  // side.

                  srnamc.SRNAMT = 'DLARHS';
                  XTYPE ??= IKL == 1 ? 'N' : 'C';
                  dlarhs(
                      PATH,
                      XTYPE,
                      'Full',
                      TRANS,
                      N,
                      N,
                      KL,
                      KU,
                      NRHS,
                      A.asMatrix(),
                      LDA,
                      XACT.asMatrix(),
                      LDB,
                      B.asMatrix(),
                      LDB,
                      ISEED,
                      INFO);
                  XTYPE = 'C';
                  dlacpy(
                      'Full', N, NRHS, B.asMatrix(), LDB, BSAV.asMatrix(), LDB);

                  if (NOFACT && ITRAN == 1) {
                    // --- Test DGBSV  ---

                    // Compute the LU factorization of the matrix
                    // and solve the system.

                    dlacpy('Full', KL + KU + 1, N, A.asMatrix(), LDA,
                        AFB(KL + 1).asMatrix(), LDAFB);
                    dlacpy(
                        'Full', N, NRHS, B.asMatrix(), LDB, X.asMatrix(), LDB);

                    srnamc.SRNAMT = 'DGBSV';
                    dgbsv(N, KL, KU, NRHS, AFB.asMatrix(), LDAFB, IWORK,
                        X.asMatrix(), LDB, INFO);

                    // Check error code from DGBSV .
                    test.expect(INFO.value, IZERO);
                    if (INFO.value != IZERO) {
                      alaerh(PATH, 'DGBSV ', INFO.value, IZERO, ' ', N, N, KL,
                          KU, NRHS, IMAT, NFAIL, NERRS, NOUT);
                    }

                    // Reconstruct matrix from factors and
                    // compute residual.

                    dgbt01(N, N, KL, KU, A.asMatrix(), LDA, AFB.asMatrix(),
                        LDAFB, IWORK, WORK, RESULT(1));
                    final int NT;
                    if (IZERO == 0) {
                      // Compute residual of the computed
                      // solution.

                      dlacpy('Full', N, NRHS, B.asMatrix(), LDB,
                          WORK.asMatrix(), LDB);
                      dgbt02(
                          'No transpose',
                          N,
                          N,
                          KL,
                          KU,
                          NRHS,
                          A.asMatrix(),
                          LDA,
                          X.asMatrix(),
                          LDB,
                          WORK.asMatrix(),
                          LDB,
                          RWORK,
                          RESULT(2));

                      // Check solution from generated exact
                      // solution.

                      dget04(N, NRHS, X.asMatrix(), LDB, XACT.asMatrix(), LDB,
                          RCONDC, RESULT(3));
                      NT = 3;
                    } else {
                      NT = 1;
                    }

                    // Print information about the tests that did
                    // not pass the threshold.

                    for (var K = 1; K <= NT; K++) {
                      final reason =
                          ' DGBSV , N=${N.i5}, KL=${KL.i5}, KU=${KU.i5}, type ${IMAT.i1}, test(${K.i1})=${RESULT[K].g12_5}';
                      test.expect(RESULT[K], lessThan(THRESH), reason: reason);
                      if (RESULT[K] >= THRESH) {
                        if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                        NOUT.println(reason);
                        NFAIL++;
                      }
                    }
                    NRUN += NT;
                  }

                  // --- Test DGBSVX ---

                  if (!PREFAC) {
                    dlaset('Full', 2 * KL + KU + 1, N, ZERO, ZERO,
                        AFB.asMatrix(), LDAFB);
                  }
                  dlaset('Full', N, NRHS, ZERO, ZERO, X.asMatrix(), LDB);
                  if (IEQUED > 1 && N > 0) {
                    // Equilibrate the matrix if FACT = 'F' and
                    // EQUED = 'R', 'C', or 'B'.

                    dlaqgb(N, N, KL, KU, A.asMatrix(), LDA, S, S(N + 1),
                        ROWCND.value, COLCND.value, AMAX.value, EQUED);
                  }

                  // Solve the system and compute the condition
                  // number and error bounds using DGBSVX.

                  final RCOND = Box(ZERO);
                  srnamc.SRNAMT = 'DGBSVX';
                  dgbsvx(
                      FACT,
                      TRANS,
                      N,
                      KL,
                      KU,
                      NRHS,
                      A.asMatrix(),
                      LDA,
                      AFB.asMatrix(),
                      LDAFB,
                      IWORK,
                      EQUED,
                      S,
                      S(N + 1),
                      B.asMatrix(),
                      LDB,
                      X.asMatrix(),
                      LDB,
                      RCOND,
                      RWORK,
                      RWORK(NRHS + 1),
                      WORK,
                      IWORK(N + 1),
                      INFO);

                  // Check the error code from DGBSVX.
                  test.expect(INFO.value, IZERO);
                  if (INFO.value != IZERO) {
                    alaerh(PATH, 'DGBSVX', INFO.value, IZERO, FACT + TRANS, N,
                        N, KL, KU, NRHS, IMAT, NFAIL, NERRS, NOUT);
                  }

                  // Compare WORK(1) from DGBSVX with the computed
                  // reciprocal pivot growth factor RPVGRW

                  double RPVGRW;
                  if (INFO.value != 0 && INFO.value <= N) {
                    var ANRMPV = ZERO;
                    for (var J = 1; J <= INFO.value; J++) {
                      for (var I = max(KU + 2 - J, 1);
                          I <= min(N + KU + 1 - J, KL + KU + 1);
                          I++) {
                        ANRMPV = max(ANRMPV, A[I + (J - 1) * LDA].abs());
                      }
                    }
                    RPVGRW = dlantb(
                        'M',
                        'U',
                        'N',
                        INFO.value,
                        min(INFO.value - 1, KL + KU),
                        AFB(max(1, KL + KU + 2 - INFO.value)).asMatrix(),
                        LDAFB,
                        WORK);
                    if (RPVGRW == ZERO) {
                      RPVGRW = ONE;
                    } else {
                      RPVGRW = ANRMPV / RPVGRW;
                    }
                  } else {
                    RPVGRW = dlantb(
                        'M', 'U', 'N', N, KL + KU, AFB.asMatrix(), LDAFB, WORK);
                    if (RPVGRW == ZERO) {
                      RPVGRW = ONE;
                    } else {
                      RPVGRW = dlangb('M', N, KL, KU, A.asMatrix(), LDA, WORK) /
                          RPVGRW;
                    }
                  }
                  RESULT[7] = (RPVGRW - WORK[1]).abs() /
                      max(WORK[1], RPVGRW) /
                      dlamch('E');

                  final int K1;
                  if (!PREFAC) {
                    // Reconstruct matrix from factors and
                    // compute residual.

                    dgbt01(N, N, KL, KU, A.asMatrix(), LDA, AFB.asMatrix(),
                        LDAFB, IWORK, WORK, RESULT(1));
                    K1 = 1;
                  } else {
                    K1 = 2;
                  }

                  final bool TRFCON;
                  if (INFO.value == 0) {
                    TRFCON = false;

                    // Compute residual of the computed solution.

                    dlacpy('Full', N, NRHS, BSAV.asMatrix(), LDB,
                        WORK.asMatrix(), LDB);
                    dgbt02(
                        TRANS,
                        N,
                        N,
                        KL,
                        KU,
                        NRHS,
                        ASAV.asMatrix(),
                        LDA,
                        X.asMatrix(),
                        LDB,
                        WORK.asMatrix(),
                        LDB,
                        RWORK(2 * NRHS + 1),
                        RESULT(2));

                    // Check solution from generated exact
                    // solution.

                    if (NOFACT || (PREFAC && lsame(EQUED.value, 'N'))) {
                      dget04(N, NRHS, X.asMatrix(), LDB, XACT.asMatrix(), LDB,
                          RCONDC, RESULT(3));
                    } else {
                      final ROLDC = ITRAN == 1 ? ROLDO : ROLDI;
                      dget04(N, NRHS, X.asMatrix(), LDB, XACT.asMatrix(), LDB,
                          ROLDC, RESULT(3));
                    }

                    // Check the error bounds from iterative
                    // refinement.

                    dgbt05(
                        TRANS,
                        N,
                        KL,
                        KU,
                        NRHS,
                        ASAV.asMatrix(),
                        LDA,
                        B.asMatrix(),
                        LDB,
                        X.asMatrix(),
                        LDB,
                        XACT.asMatrix(),
                        LDB,
                        RWORK,
                        RWORK(NRHS + 1),
                        RESULT(4));
                  } else {
                    TRFCON = true;
                  }

                  // Compare RCOND from DGBSVX with the computed
                  // value in RCONDC.

                  RESULT[6] = dget06(RCOND.value, RCONDC);

                  // Print information about the tests that did
                  // not pass the threshold.

                  final tests = !TRFCON
                      ? [for (var K = K1; K <= NTESTS; K++) K]
                      : [if (!PREFAC) 1, 6, 7];
                  for (var K in tests) {
                    final reason = PREFAC
                        ? ' DGBSVX( \'${FACT.a1}\',\'${TRANS.a1}\',${N.i5},${KL.i5},${KU.i5},...), EQUED=\'${EQUED.value.a1}\', type ${IMAT.i1}, test(${K.i1})=${RESULT[K].g12_5}'
                        : ' DGBSVX( \'${FACT.a1}\',\'${TRANS.a1}\',${N.i5},${KL.i5},${KU.i5},...), type ${IMAT.i1}, test(${K.i1})=${RESULT[K].g12_5}';
                    test.expect(RESULT[K], lessThan(THRESH), reason: reason);
                    if (RESULT[K] >= THRESH) {
                      if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                      NOUT.println(reason);
                      NFAIL++;
                      NRUN++;
                    }
                  }
                }
              }
            }
          }, skip: skip);
        }
      }
    }
  }

  // Print a summary of the results.
  alasvm(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
