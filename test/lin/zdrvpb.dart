import 'dart:math';

import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlanhb.dart';
import 'package:lapack/src/zlaqhb.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zpbequ.dart';
import 'package:lapack/src/zpbsv.dart';
import 'package:lapack/src/zpbsvx.dart';
import 'package:lapack/src/zpbtrf.dart';
import 'package:lapack/src/zpbtrs.dart';

import '../matgen/zlatms.dart';
import 'aladhd.dart';
import 'alaerh.dart';
import 'alasvm.dart';
import 'common.dart';
import 'dget06.dart';
import 'xlaenv.dart';
import 'zerrvxx.dart';
import 'zget04.dart';
import 'zlaipd.dart';
import 'zlarhs.dart';
import 'zlatb4.dart';
import 'zpbt01.dart';
import 'zpbt02.dart';
import 'zpbt05.dart';

void zdrvpb(
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
  const NTYPES = 8, NTESTS = 6;
  const NBW = 4;
  final ISEED = Array<int>(4), KDVAL = Array<int>(NBW);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const FACTS = ['F', 'N', 'E'], EQUEDS = ['N', 'Y'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Zomplex precision'[0]}PB';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) zerrvx(PATH, NOUT);
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

        var IZERO = 0, I1 = 0, I2 = 0;
        for (var IMAT = 1; IMAT <= NIMAT; IMAT++) {
          // Do the tests only if DOTYPE( IMAT ) is true.

          if (!DOTYPE[IMAT]) continue;

          // Skip types 2, 3, or 4 if the matrix size is too small.

          final ZEROT = IMAT >= 2 && IMAT <= 4;
          if (ZEROT && N < IMAT - 1) continue;

          if (!ZEROT || !DOTYPE[1]) {
            // Set up parameters with ZLATB4 and generate a test
            // matrix with ZLATMS.

            final (:TYPE, KL: _, KU: _, :ANORM, :MODE, :CNDNUM, :DIST) =
                zlatb4(PATH, IMAT, N, N);

            srnamc.SRNAMT = 'ZLATMS';
            zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KD, KD,
                PACKIT, A(KOFF).asMatrix(), LDAB, WORK, INFO);

            // Check error code from ZLATMS.

            if (INFO.value != 0) {
              alaerh(PATH, 'ZLATMS', INFO.value, 0, UPLO, N, N, -1, -1, -1,
                  IMAT, NFAIL, NERRS, NOUT);
              continue;
            }
          } else if (IZERO > 0) {
            // Use the same matrix for types 3 and 4 as for type
            // 2 by copying back the zeroed out column,

            var IW = 2 * LDA + 1;
            if (IUPLO == 1) {
              final IOFF = (IZERO - 1) * LDAB + KD + 1;
              zcopy(IZERO - I1, WORK(IW), 1, A(IOFF - IZERO + I1), 1);
              IW += IZERO - I1;
              zcopy(I2 - IZERO + 1, WORK(IW), 1, A(IOFF), max(LDAB - 1, 1));
            } else {
              var IOFF = (I1 - 1) * LDAB + 1;
              zcopy(IZERO - I1, WORK(IW), 1, A(IOFF + IZERO - I1),
                  max(LDAB - 1, 1));
              IOFF = (IZERO - 1) * LDAB + 1;
              IW += IZERO - I1;
              zcopy(I2 - IZERO + 1, WORK(IW), 1, A(IOFF), 1);
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
              WORK[IW + I] = Complex.zero;
            }
            IW++;
            final I1 = max(IZERO - KD, 1);
            final I2 = min(IZERO + KD, N);

            if (IUPLO == 1) {
              final IOFF = (IZERO - 1) * LDAB + KD + 1;
              zswap(IZERO - I1, A(IOFF - IZERO + I1), 1, WORK(IW), 1);
              IW += IZERO - I1;
              zswap(I2 - IZERO + 1, A(IOFF), max(LDAB - 1, 1), WORK(IW), 1);
            } else {
              var IOFF = (I1 - 1) * LDAB + 1;
              zswap(IZERO - I1, A(IOFF + IZERO - I1), max(LDAB - 1, 1),
                  WORK(IW), 1);
              IOFF = (IZERO - 1) * LDAB + 1;
              IW += IZERO - I1;
              zswap(I2 - IZERO + 1, A(IOFF), 1, WORK(IW), 1);
            }
          }

          // Set the imaginary part of the diagonals.

          if (IUPLO == 1) {
            zlaipd(N, A(KD + 1), LDAB, 0);
          } else {
            zlaipd(N, A(1), LDAB, 0);
          }

          // Save a copy of the matrix A in ASAV.

          zlacpy('Full', KD + 1, N, A.asMatrix(), LDAB, ASAV.asMatrix(), LDAB);

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
                // Compute the condition number for comparison
                // with the value returned by ZPBSVX (FACT =
                // 'N' reuses the condition number from the
                // previous iteration with FACT = 'F').

                zlacpy('Full', KD + 1, N, ASAV.asMatrix(), LDAB,
                    AFAC.asMatrix(), LDAB);
                if (EQUIL || IEQUED > 1) {
                  // Compute row and column scale factors to
                  // equilibrate the matrix A.

                  zpbequ(
                      UPLO, N, KD, AFAC.asMatrix(), LDAB, S, SCOND, AMAX, INFO);
                  if (INFO.value == 0 && N > 0) {
                    if (IEQUED > 1) SCOND.value = ZERO;

                    // Equilibrate the matrix.

                    zlaqhb(UPLO, N, KD, AFAC.asMatrix(), LDAB, S, SCOND.value,
                        AMAX.value, EQUED);
                  }
                }

                // Save the condition number of the
                // non-equilibrated system for use in ZGET04.

                if (EQUIL) ROLDC = RCONDC;

                // Compute the 1-norm of A.

                final ANORM =
                    zlanhb('1', UPLO, N, KD, AFAC.asMatrix(), LDAB, RWORK);

                // Factor the matrix A.

                zpbtrf(UPLO, N, KD, AFAC.asMatrix(), LDAB, INFO);

                // Form the inverse of A.

                zlaset(
                    'Full', N, N, Complex.zero, Complex.one, A.asMatrix(), LDA);
                srnamc.SRNAMT = 'ZPBTRS';
                zpbtrs(UPLO, N, KD, N, AFAC.asMatrix(), LDAB, A.asMatrix(), LDA,
                    INFO);

                // Compute the 1-norm condition number of A.

                final AINVNM = zlange('1', N, N, A.asMatrix(), LDA, RWORK);
                if (ANORM <= ZERO || AINVNM <= ZERO) {
                  RCONDC = ONE;
                } else {
                  RCONDC = (ONE / ANORM) / AINVNM;
                }
              }

              // Restore the matrix A.

              zlacpy(
                  'Full', KD + 1, N, ASAV.asMatrix(), LDAB, A.asMatrix(), LDAB);

              // Form an exact solution and set the right hand
              // side.

              srnamc.SRNAMT = 'ZLARHS';
              zlarhs(PATH, XTYPE, UPLO, ' ', N, N, KD, KD, NRHS, A.asMatrix(),
                  LDAB, XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
              XTYPE = 'C';
              zlacpy('Full', N, NRHS, B.asMatrix(), LDA, BSAV.asMatrix(), LDA);

              if (NOFACT) {
                // --- Test ZPBSV  ---

                // Compute the L*L' or U'*U factorization of the
                // matrix and solve the system.

                zlacpy('Full', KD + 1, N, A.asMatrix(), LDAB, AFAC.asMatrix(),
                    LDAB);
                zlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

                srnamc.SRNAMT = 'ZPBSV ';
                zpbsv(UPLO, N, KD, NRHS, AFAC.asMatrix(), LDAB, X.asMatrix(),
                    LDA, INFO);

                // Check error code from ZPBSV .

                if (INFO.value != IZERO) {
                  alaerh(PATH, 'ZPBSV ', INFO.value, IZERO, UPLO, N, N, KD, KD,
                      NRHS, IMAT, NFAIL, NERRS, NOUT);
                } else if (INFO.value != 0) {
                  //
                } else {
                  // Reconstruct matrix from factors and compute
                  // residual.

                  zpbt01(UPLO, N, KD, A.asMatrix(), LDAB, AFAC.asMatrix(), LDAB,
                      RWORK, RESULT(1));

                  // Compute residual of the computed solution.

                  zlacpy(
                      'Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
                  zpbt02(UPLO, N, KD, NRHS, A.asMatrix(), LDAB, X.asMatrix(),
                      LDA, WORK.asMatrix(), LDA, RWORK, RESULT(2));

                  // Check solution from generated exact solution.

                  zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                      RCONDC, RESULT(3));
                  const NT = 3;

                  // Print information about the tests that did
                  // not pass the threshold.

                  for (var K = 1; K <= NT; K++) {
                    if (RESULT[K] >= THRESH) {
                      if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                      NOUT.println(
                          ' ZPBSV , UPLO=\'${UPLO.a1}\', N =${N.i5}, KD =${KD.i5}, type ${IMAT.i1}, test(${K.i1})=${RESULT[K].g12_5}');
                      NFAIL++;
                    }
                  }
                  NRUN += NT;
                }
              }

              // --- Test ZPBSVX ---

              if (!PREFAC) {
                zlaset('Full', KD + 1, N, Complex.zero, Complex.zero,
                    AFAC.asMatrix(), LDAB);
              }
              zlaset('Full', N, NRHS, Complex.zero, Complex.zero, X.asMatrix(),
                  LDA);
              if (IEQUED > 1 && N > 0) {
                // Equilibrate the matrix if FACT='F' and
                // EQUED='Y'

                zlaqhb(UPLO, N, KD, A.asMatrix(), LDAB, S, SCOND.value,
                    AMAX.value, EQUED);
              }

              // Solve the system and compute the condition
              // number and error bounds using ZPBSVX.

              final RCOND = Box(ZERO);
              srnamc.SRNAMT = 'ZPBSVX';
              zpbsvx(
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
                  RWORK(2 * NRHS + 1),
                  INFO);

              // Check the error code from ZPBSVX.

              if (INFO.value != IZERO) {
                alaerh(PATH, 'ZPBSVX', INFO.value, IZERO, FACT + UPLO, N, N, KD,
                    KD, NRHS, IMAT, NFAIL, NERRS, NOUT);
                continue;
              }

              final int K1;
              if (INFO.value == 0) {
                if (!PREFAC) {
                  // Reconstruct matrix from factors and
                  // compute residual.

                  zpbt01(UPLO, N, KD, A.asMatrix(), LDAB, AFAC.asMatrix(), LDAB,
                      RWORK(2 * NRHS + 1), RESULT(1));
                  K1 = 1;
                } else {
                  K1 = 2;
                }

                // Compute residual of the computed solution.

                zlacpy('Full', N, NRHS, BSAV.asMatrix(), LDA, WORK.asMatrix(),
                    LDA);
                zpbt02(UPLO, N, KD, NRHS, ASAV.asMatrix(), LDAB, X.asMatrix(),
                    LDA, WORK.asMatrix(), LDA, RWORK(2 * NRHS + 1), RESULT(2));

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

                zpbt05(
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

              // Compare RCOND from ZPBSVX with the computed
              // value in RCONDC.

              RESULT[6] = dget06(RCOND.value, RCONDC);

              // Print information about the tests that did not
              // pass the threshold.

              for (var K = K1; K <= 6; K++) {
                if (RESULT[K] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                  if (PREFAC) {
                    NOUT.println(
                        ' ZPBSVX( \'${FACT.a1}\'${UPLO.a1}\', ${N.i5}, ${KD.i5}, ... ), EQUED=\'${EQUED.value.a1}\', type ${IMAT.i1}, test(${K.i1})=${RESULT[K].g12_5}');
                  } else {
                    NOUT.println(
                        ' ZPBSVX( \'${FACT.a1}\'${UPLO.a1}\', ${N.i5}, ${KD.i5}, ... ), type ${IMAT.i1}, test(${K.i1})=${RESULT[K].g12_5}');
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
