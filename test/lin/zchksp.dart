import 'dart:math';

import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlansp.dart';
import 'package:lapack/src/zspcon.dart';
import 'package:lapack/src/zsprfs.dart';
import 'package:lapack/src/zsptrf.dart';
import 'package:lapack/src/zsptri.dart';
import 'package:lapack/src/zsptrs.dart';

import '../matgen/zlatms.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'dget06.dart';
import 'zerrsyx.dart';
import 'zget04.dart';
import 'zlarhs.dart';
import 'zlatb4.dart';
import 'zlatsp.dart';
import 'zppt05.dart';
import 'zspt01.dart';
import 'zspt02.dart';
import 'zspt03.dart';

void zchksp(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NNS,
  final Array<int> NSVAL_,
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
  final NSVAL = NSVAL_.having();
  final A = A_.having();
  final AFAC = AFAC_.having();
  final AINV = AINV_.having();
  final B = B_.having();
  final X = X_.having();
  final XACT = XACT_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0;
  const NTYPES = 11, NTESTS = 8;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Zomplex precision'[0]}SP';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) zerrsy(PATH, NOUT);
  infoc.INFOT = 0;

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
        final PACKIT = lsame(UPLO, 'U') ? 'C' : 'R';

        var IZERO = 0, KL = 0, KU = 0;
        if (IMAT != NTYPES) {
          // Set up parameters with ZLATB4 and generate a test
          // matrix with ZLATMS.

          final (:TYPE, KL: KLTMP, KU: KUTMP, :ANORM, :MODE, :CNDNUM, :DIST) =
              zlatb4(PATH, IMAT, N, N);
          (KL, KU) = (KLTMP, KUTMP);

          srnamc.SRNAMT = 'ZLATMS';
          zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
              PACKIT, A.asMatrix(), LDA, WORK, INFO);

          // Check error code from ZLATMS.

          if (INFO.value != 0) {
            alaerh(PATH, 'ZLATMS', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
                NFAIL, NERRS, NOUT);
            continue;
          }

          // For types 3-6, zero one or more rows and columns of
          // the matrix to test that INFO.value is returned correctly.

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
              if (IUPLO == 1) {
                // Set the first IZERO rows and columns to zero.

                var IOFF = 0;
                for (var J = 1; J <= N; J++) {
                  final I2 = min(J, IZERO);
                  for (var I = 1; I <= I2; I++) {
                    A[IOFF + I] = Complex.zero;
                  }
                  IOFF += J;
                }
              } else {
                // Set the last IZERO rows and columns to zero.

                var IOFF = 0;
                for (var J = 1; J <= N; J++) {
                  final I1 = max(J, IZERO);
                  for (var I = I1; I <= N; I++) {
                    A[IOFF + I] = Complex.zero;
                  }
                  IOFF += N - J;
                }
              }
            }
          } else {
            IZERO = 0;
          }
        } else {
          // Use a special block diagonal matrix to test alternate
          // code for the 2 x 2 blocks.

          zlatsp(UPLO, N, A, ISEED);
        }

        // Compute the L*D*L' or U*D*U' factorization of the matrix.

        final NPP = N * (N + 1) ~/ 2;
        zcopy(NPP, A, 1, AFAC, 1);
        srnamc.SRNAMT = 'ZSPTRF';
        zsptrf(UPLO, N, AFAC, IWORK, INFO);

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

        // Check error code from ZSPTRF.

        if (INFO.value != K) {
          alaerh(PATH, 'ZSPTRF', INFO.value, K, UPLO, N, N, -1, -1, -1, IMAT,
              NFAIL, NERRS, NOUT);
        }
        final TRFCON = INFO.value != 0;

        // +    TEST 1
        // Reconstruct matrix from factors and compute residual.

        zspt01(UPLO, N, A, AFAC, IWORK, AINV.asMatrix(), LDA, RWORK, RESULT(1));

        // +    TEST 2
        // Form the inverse and compute the residual.

        final int NT;
        final RCONDC = Box(ZERO);
        if (!TRFCON) {
          zcopy(NPP, AFAC, 1, AINV, 1);
          srnamc.SRNAMT = 'ZSPTRI';
          zsptri(UPLO, N, AINV, IWORK, WORK, INFO);

          // Check error code from ZSPTRI.

          if (INFO.value != 0) {
            alaerh(PATH, 'ZSPTRI', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
                NFAIL, NERRS, NOUT);
          }

          zspt03(
              UPLO, N, A, AINV, WORK.asMatrix(), LDA, RWORK, RCONDC, RESULT(2));
          NT = 2;
        } else {
          NT = 1;
        }

        // Print information about the tests that did not pass
        // the threshold.

        for (var K = 1; K <= NT; K++) {
          if (RESULT[K] >= THRESH) {
            if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
            NOUT.println(
                ' UPLO = \'${UPLO.a1}\', N =${N.i5}, type ${IMAT.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}');
            NFAIL++;
          }
        }
        NRUN += NT;

        // Do only the condition estimate if INFO.value is not 0.

        if (TRFCON) {
          RCONDC.value = ZERO;
        } else {
          for (var IRHS = 1; IRHS <= NNS; IRHS++) {
            final NRHS = NSVAL[IRHS];

            // +    TEST 3
            // Solve and compute residual for  A * X = B.

            srnamc.SRNAMT = 'ZLARHS';
            zlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A.asMatrix(),
                LDA, XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
            zlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

            srnamc.SRNAMT = 'ZSPTRS';
            zsptrs(UPLO, N, NRHS, AFAC, IWORK, X.asMatrix(), LDA, INFO);

            // Check error code from ZSPTRS.

            if (INFO.value != 0) {
              alaerh(PATH, 'ZSPTRS', INFO.value, 0, UPLO, N, N, -1, -1, NRHS,
                  IMAT, NFAIL, NERRS, NOUT);
            }

            zlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
            zspt02(UPLO, N, NRHS, A, X.asMatrix(), LDA, WORK.asMatrix(), LDA,
                RWORK, RESULT(3));

            // +    TEST 4
            // Check solution from generated exact solution.

            zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                RCONDC.value, RESULT(4));

            // +    TESTS 5, 6, and 7
            // Use iterative refinement to improve the solution.

            srnamc.SRNAMT = 'ZSPRFS';
            zsprfs(
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
                RWORK,
                RWORK(NRHS + 1),
                WORK,
                RWORK(2 * NRHS + 1),
                INFO);

            // Check error code from ZSPRFS.

            if (INFO.value != 0) {
              alaerh(PATH, 'ZSPRFS', INFO.value, 0, UPLO, N, N, -1, -1, NRHS,
                  IMAT, NFAIL, NERRS, NOUT);
            }

            zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                RCONDC.value, RESULT(5));
            zppt05(UPLO, N, NRHS, A, B.asMatrix(), LDA, X.asMatrix(), LDA,
                XACT.asMatrix(), LDA, RWORK, RWORK(NRHS + 1), RESULT(6));

            // Print information about the tests that did not pass
            // the threshold.

            for (var K = 3; K <= 7; K++) {
              if (RESULT[K] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                NOUT.println(
                    ' UPLO = \'${UPLO.a1}\', N =${N.i5}, NRHS=${NRHS.i3}, type ${IMAT.i2}, test(${K.i2}) =${RESULT[K].g12_5}');
                NFAIL++;
              }
            }
            NRUN += 5;
          }
        }

        // +    TEST 8
        // Get an estimate of RCOND = 1/CNDNUM.

        final ANORM = zlansp('1', UPLO, N, A, RWORK);
        srnamc.SRNAMT = 'ZSPCON';
        final RCOND = Box(ZERO);
        zspcon(UPLO, N, AFAC, IWORK, ANORM, RCOND, WORK, INFO);

        // Check error code from ZSPCON.

        if (INFO.value != 0) {
          alaerh(PATH, 'ZSPCON', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
              NFAIL, NERRS, NOUT);
        }

        RESULT[8] = dget06(RCOND.value, RCONDC.value);

        // Print the test ratio if it is >= THRESH.

        if (RESULT[8] >= THRESH) {
          if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
          NOUT.println(
              ' UPLO = \'${UPLO.a1}\', N =${N.i5}, type ${IMAT.i2}, test ${8.i2}, ratio =${RESULT[8].g12_5}');
          NFAIL++;
        }
        NRUN++;
      }
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
