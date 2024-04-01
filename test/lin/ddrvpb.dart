import 'dart:math';

import 'package:lapack/lapack.dart';

import '../matgen/dlatms.dart';
import '../test_driver.dart';
import 'aladhd.dart';
import 'alaerh.dart';
import 'alasvm.dart';
import 'common.dart';
import 'derrvx.dart';
import 'dget04.dart';
import 'dget06.dart';
import 'dlarhs.dart';
import 'dlatb4.dart';
import 'dpbt01.dart';
import 'dpbt02.dart';
import 'dpbt05.dart';
import 'xlaenv.dart';

void ddrvpb(
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
  final TestDriver test,
) {
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

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0, ZERO = 0.0;
  const NTYPES = 8, NTESTS = 6;
  const NBW = 4;
  final ISEED = Array<int>(4), KDVAL = Array<int>(NBW);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const FACTS = ['F', 'N', 'E'];
  const EQUEDS = ['N', 'Y'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}PB';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) derrvx(PATH, NOUT, test);
  infoc.INFOT = 0;
  KDVAL[1] = 0;

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

    // Set limits on the number of loop iterations.

    final NKD = max(1, min(N, 4));
    final NIMAT = N == 0 ? 1 : NTYPES;

    KDVAL[2] = N + (N + 1) ~/ 4;
    KDVAL[3] = (3 * N - 1) ~/ 4;
    KDVAL[4] = (N + 1) ~/ 4;

    for (var IKD = 1; IKD <= NKD; IKD++) {
      // Do for KD = 0, (5*N+1)/4, (3N-1)/4, and (N+1)/4. This order
      // makes it easier to skip redundant values for small values
      // of N.

      final KD = KDVAL[IKD];
      final LDAB = KD + 1;

      // Do first for UPLO = 'U', then for UPLO = 'L'

      for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
        final (UPLO, PACKIT, KOFF) =
            IUPLO == 1 ? ('U', 'Q', max(1, KD + 2 - N)) : ('L', 'B', 1);

        int IZERO = 0, I1 = 0, I2 = 0;
        for (var IMAT = 1; IMAT <= NIMAT; IMAT++) {
          // Do the tests only if DOTYPE( IMAT ) is true.

          if (!DOTYPE[IMAT]) continue;

          // Skip types 2, 3, or 4 if the matrix size is too small.

          final ZEROT = IMAT >= 2 && IMAT <= 4;
          if (ZEROT && N < IMAT - 1) continue;

          if (!ZEROT || !DOTYPE[1]) {
            // Set up parameters with DLATB4 and generate a test
            // matrix with DLATMS.
            final (:TYPE, KL: _, KU: _, :ANORM, :MODE, :COND, :DIST) =
                dlatb4(PATH, IMAT, N, N);

            srnamc.SRNAMT = 'DLATMS';
            dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, COND, ANORM, KD, KD,
                PACKIT, A(KOFF).asMatrix(), LDAB, WORK, INFO);

            // Check error code from DLATMS.

            if (INFO.value != 0) {
              alaerh(PATH, 'DLATMS', INFO.value, 0, UPLO, N, N, -1, -1, -1,
                  IMAT, NFAIL, NERRS, NOUT);
              continue;
            }
          } else if (IZERO > 0) {
            // Use the same matrix for types 3 and 4 as for type
            // 2 by copying back the zeroed out column,

            var IW = 2 * LDA + 1;
            if (IUPLO == 1) {
              final IOFF = (IZERO - 1) * LDAB + KD + 1;
              dcopy(IZERO - I1, WORK(IW), 1, A(IOFF - IZERO + I1), 1);
              IW += IZERO - I1;
              dcopy(I2 - IZERO + 1, WORK(IW), 1, A(IOFF), max(LDAB - 1, 1));
            } else {
              var IOFF = (I1 - 1) * LDAB + 1;
              dcopy(IZERO - I1, WORK(IW), 1, A(IOFF + IZERO - I1),
                  max(LDAB - 1, 1));
              IOFF = (IZERO - 1) * LDAB + 1;
              IW += IZERO - I1;
              dcopy(I2 - IZERO + 1, WORK(IW), 1, A(IOFF), 1);
            }
          }

          // For types 2-4, zero one row and column of the matrix
          // to test that INFO is returned correctly.

          IZERO = 0;
          if (ZEROT) {
            if (IMAT == 2) {
              IZERO = 1;
            } else if (IMAT == 3) {
              IZERO = N;
            } else {
              IZERO = N ~/ 2 + 1;
            }

            // Save the zeroed out row and column in WORK(*,3)

            var IW = 2 * LDA;
            for (var I = 1; I <= min(2 * KD + 1, N); I++) {
              WORK[IW + I] = ZERO;
            }
            IW++;
            I1 = max(IZERO - KD, 1);
            I2 = min(IZERO + KD, N);

            if (IUPLO == 1) {
              final IOFF = (IZERO - 1) * LDAB + KD + 1;
              dswap(IZERO - I1, A(IOFF - IZERO + I1), 1, WORK(IW), 1);
              IW += IZERO - I1;
              dswap(I2 - IZERO + 1, A(IOFF), max(LDAB - 1, 1), WORK(IW), 1);
            } else {
              var IOFF = (I1 - 1) * LDAB + 1;
              dswap(IZERO - I1, A(IOFF + IZERO - I1), max(LDAB - 1, 1),
                  WORK(IW), 1);
              IOFF = (IZERO - 1) * LDAB + 1;
              IW += IZERO - I1;
              dswap(I2 - IZERO + 1, A(IOFF), 1, WORK(IW), 1);
            }
          }

          // Save a copy of the matrix A in ASAV.

          dlacpy('Full', KD + 1, N, A.asMatrix(), LDAB, ASAV.asMatrix(), LDAB);

          for (var IEQUED = 1; IEQUED <= 2; IEQUED++) {
            final EQUED = Box(EQUEDS[IEQUED - 1]);
            final NFACT = IEQUED == 1 ? 3 : 1;

            var RCONDC = ZERO, ROLDC = ZERO;
            final SCOND = Box(0.0), AMAX = Box(0.0);
            for (var IFACT = 1; IFACT <= NFACT; IFACT++) {
              final FACT = FACTS[IFACT - 1];
              final PREFAC = lsame(FACT, 'F');
              final NOFACT = lsame(FACT, 'N');
              final EQUIL = lsame(FACT, 'E');

              if (ZEROT) {
                if (PREFAC) continue;
                RCONDC = ZERO;
              } else if (!lsame(FACT, 'N')) {
                // Compute the condition number for comparison
                // with the value returned by DPBSVX (FACT =
                // 'N' reuses the condition number from the
                // previous iteration with FACT = 'F').

                dlacpy('Full', KD + 1, N, ASAV.asMatrix(), LDAB,
                    AFAC.asMatrix(), LDAB);
                if (EQUIL || IEQUED > 1) {
                  // Compute row and column scale factors to
                  // equilibrate the matrix A.

                  dpbequ(
                      UPLO, N, KD, AFAC.asMatrix(), LDAB, S, SCOND, AMAX, INFO);
                  if (INFO.value == 0 && N > 0) {
                    if (IEQUED > 1) SCOND.value = ZERO;

                    // Equilibrate the matrix.

                    dlaqsb(UPLO, N, KD, AFAC.asMatrix(), LDAB, S, SCOND.value,
                        AMAX.value, EQUED);
                  }
                }

                // Save the condition number of the
                // non-equilibrated system for use in DGET04.

                if (EQUIL) ROLDC = RCONDC;

                // Compute the 1-norm of A.

                final ANORM =
                    dlansb('1', UPLO, N, KD, AFAC.asMatrix(), LDAB, RWORK);

                // Factor the matrix A.

                dpbtrf(UPLO, N, KD, AFAC.asMatrix(), LDAB, INFO);

                // Form the inverse of A.

                dlaset('Full', N, N, ZERO, ONE, A.asMatrix(), LDA);
                srnamc.SRNAMT = 'DPBTRS';
                dpbtrs(UPLO, N, KD, N, AFAC.asMatrix(), LDAB, A.asMatrix(), LDA,
                    INFO);

                // Compute the 1-norm condition number of A.

                final AINVNM = dlange('1', N, N, A.asMatrix(), LDA, RWORK);
                if (ANORM <= ZERO || AINVNM <= ZERO) {
                  RCONDC = ONE;
                } else {
                  RCONDC = (ONE / ANORM) / AINVNM;
                }
              }

              // Restore the matrix A.

              dlacpy(
                  'Full', KD + 1, N, ASAV.asMatrix(), LDAB, A.asMatrix(), LDAB);

              // Form an exact solution and set the right hand
              // side.

              srnamc.SRNAMT = 'DLARHS';
              dlarhs(PATH, XTYPE, UPLO, ' ', N, N, KD, KD, NRHS, A.asMatrix(),
                  LDAB, XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
              XTYPE = 'C';
              dlacpy('Full', N, NRHS, B.asMatrix(), LDA, BSAV.asMatrix(), LDA);

              if (NOFACT) {
                // --- Test DPBSV  ---

                // Compute the L*L' or U'*U factorization of the
                // matrix and solve the system.

                dlacpy('Full', KD + 1, N, A.asMatrix(), LDAB, AFAC.asMatrix(),
                    LDAB);
                dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

                srnamc.SRNAMT = 'DPBSV ';
                dpbsv(UPLO, N, KD, NRHS, AFAC.asMatrix(), LDAB, X.asMatrix(),
                    LDA, INFO);

                // Check error code from DPBSV .

                if (INFO.value != IZERO) {
                  alaerh(PATH, 'DPBSV ', INFO.value, IZERO, UPLO, N, N, KD, KD,
                      NRHS, IMAT, NFAIL, NERRS, NOUT);
                } else if (INFO.value == 0) {
                  // Reconstruct matrix from factors and compute
                  // residual.

                  dpbt01(UPLO, N, KD, A.asMatrix(), LDAB, AFAC.asMatrix(), LDAB,
                      RWORK, RESULT(1));

                  // Compute residual of the computed solution.

                  dlacpy(
                      'Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
                  dpbt02(UPLO, N, KD, NRHS, A.asMatrix(), LDAB, X.asMatrix(),
                      LDA, WORK.asMatrix(), LDA, RWORK, RESULT(2));

                  // Check solution from generated exact solution.

                  dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                      RCONDC, RESULT(3));
                  const NT = 3;

                  // Print information about the tests that did
                  // not pass the threshold.

                  for (var K = 1; K <= NT; K++) {
                    if (RESULT[K] >= THRESH) {
                      if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                      NOUT.println(
                          ' DPBSV, UPLO=\'${UPLO.a1}\', N =${N.i5}, KD =${KD.i5}, type ${IMAT.i1}, test(${K.i1})=${RESULT[K].g12_5}');
                      NFAIL++;
                    }
                  }
                  NRUN += NT;
                }
              }

              // --- Test DPBSVX ---

              if (!PREFAC) {
                dlaset('Full', KD + 1, N, ZERO, ZERO, AFAC.asMatrix(), LDAB);
              }
              dlaset('Full', N, NRHS, ZERO, ZERO, X.asMatrix(), LDA);
              if (IEQUED > 1 && N > 0) {
                // Equilibrate the matrix if FACT='F' and
                // EQUED='Y'

                dlaqsb(UPLO, N, KD, A.asMatrix(), LDAB, S, SCOND.value,
                    AMAX.value, EQUED);
              }

              // Solve the system and compute the condition
              // number and error bounds using DPBSVX.

              final RCOND = Box(0.0);
              srnamc.SRNAMT = 'DPBSVX';
              dpbsvx(
                  FACT,
                  UPLO,
                  N,
                  KD,
                  NRHS,
                  A.asMatrix(),
                  LDAB,
                  AFAC.asMatrix(),
                  LDAB,
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

              // Check the error code from DPBSVX.

              if (INFO.value != IZERO) {
                alaerh(PATH, 'DPBSVX', INFO.value, IZERO, FACT + UPLO, N, N, KD,
                    KD, NRHS, IMAT, NFAIL, NERRS, NOUT);
                continue;
              }

              final int K1;
              if (INFO.value == 0) {
                if (!PREFAC) {
                  // Reconstruct matrix from factors and
                  // compute residual.

                  dpbt01(UPLO, N, KD, A.asMatrix(), LDAB, AFAC.asMatrix(), LDAB,
                      RWORK(2 * NRHS + 1), RESULT(1));
                  K1 = 1;
                } else {
                  K1 = 2;
                }

                // Compute residual of the computed solution.

                dlacpy('Full', N, NRHS, BSAV.asMatrix(), LDA, WORK.asMatrix(),
                    LDA);
                dpbt02(UPLO, N, KD, NRHS, ASAV.asMatrix(), LDAB, X.asMatrix(),
                    LDA, WORK.asMatrix(), LDA, RWORK(2 * NRHS + 1), RESULT(2));

                // Check solution from generated exact solution.

                if (NOFACT || (PREFAC && lsame(EQUED.value, 'N'))) {
                  dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                      RCONDC, RESULT(3));
                } else {
                  dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                      ROLDC, RESULT(3));
                }

                // Check the error bounds from iterative
                // refinement.

                dpbt05(
                    UPLO,
                    N,
                    KD,
                    NRHS,
                    ASAV.asMatrix(),
                    LDAB,
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

              // Compare RCOND from DPBSVX with the computed
              // value in RCONDC.

              RESULT[6] = dget06(RCOND.value, RCONDC);

              // Print information about the tests that did not
              // pass the threshold.

              for (var K = K1; K <= 6; K++) {
                if (RESULT[K] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                  if (PREFAC) {
                    NOUT.println(
                        ' DPBSVX( \'${FACT.a1}\'${UPLO.a1}\', ${N.i5}, ${KD.i5}, ... ), EQUED=\'${EQUED.value.a1}\', type ${IMAT.i1}, test(${K.i1})=${RESULT[K].g12_5}');
                  } else {
                    NOUT.println(
                        ' DPBSVX( \'${FACT.a1}\'${UPLO.a1}\', ${N.i5}, ${KD.i5}, ... ), type ${IMAT.i1}, test(${K.i1})=${RESULT[K].g12_5}');
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
