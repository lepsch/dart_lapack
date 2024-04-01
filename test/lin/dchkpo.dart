import 'dart:math';

import 'package:lapack/lapack.dart';

import '../matgen/dlatms.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'derrpo.dart';
import 'dget04.dart';
import 'dget06.dart';
import 'dlarhs.dart';
import 'dlatb4.dart';
import 'dpot01.dart';
import 'dpot02.dart';
import 'dpot03.dart';
import 'dpot05.dart';
import 'xlaenv.dart';

void dchkpo(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NNB,
  final Array<int> NBVAL_,
  final int NNS,
  final Array<int> NSVAL_,
  final double THRESH,
  final bool TSTERR,
  final int NMAX,
  final Array<double> A_,
  final Array<double> AFAC_,
  final Array<double> AINV_,
  final Array<double> B_,
  final Array<double> X_,
  final Array<double> XACT_,
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
  final NBVAL = NBVAL_.having();
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
  const NTYPES = 9;
  const NTESTS = 8;
  int NFAIL, NRUN;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'];
  final INFO = Box(0), NERRS = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}PO';
  NRUN = 0;
  NFAIL = 0;
  NERRS.value = 0;
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) derrpo(PATH, NOUT);
  infoc.INFOT = 0;
  xlaenv(2, 2);

  // Do for each value of N in NVAL

  for (var IN = 1; IN <= NN; IN++) {
    final N = NVAL[IN];
    final LDA = max(N, 1).toInt();
    final XTYPE = 'N';
    final NIMAT = N <= 0 ? 1 : NTYPES;

    var IZERO = 0;
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
        // test that INFO is returned correctly.

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

        // Do for each value of NB in NBVAL

        for (var INB = 1; INB <= NNB; INB++) {
          final NB = NBVAL[INB];
          xlaenv(1, NB);

          // Compute the L*L' or U'*U factorization of the matrix.

          dlacpy(UPLO, N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);
          srnamc.SRNAMT = 'DPOTRF';
          dpotrf(UPLO, N, AFAC.asMatrix(), LDA, INFO);

          // Check error code from DPOTRF.

          if (INFO.value != IZERO) {
            alaerh(PATH, 'DPOTRF', INFO.value, IZERO, UPLO, N, N, -1, -1, NB,
                IMAT, NFAIL, NERRS, NOUT);
            continue;
          }

          // Skip the tests if INFO is not 0.

          if (INFO.value != 0) continue;

          // +    TEST 1
          // Reconstruct matrix from factors and compute residual.

          dlacpy(UPLO, N, N, AFAC.asMatrix(), LDA, AINV.asMatrix(), LDA);
          dpot01(UPLO, N, A.asMatrix(), LDA, AINV.asMatrix(), LDA, RWORK,
              RESULT(1));

          // +    TEST 2
          // Form the inverse and compute the residual.

          dlacpy(UPLO, N, N, AFAC.asMatrix(), LDA, AINV.asMatrix(), LDA);
          srnamc.SRNAMT = 'DPOTRI';
          dpotri(UPLO, N, AINV.asMatrix(), LDA, INFO);

          // Check error code from DPOTRI.

          if (INFO.value != 0) {
            alaerh(PATH, 'DPOTRI', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
                NFAIL, NERRS, NOUT);
          }

          final RCONDC = Box(0.0);
          dpot03(UPLO, N, A.asMatrix(), LDA, AINV.asMatrix(), LDA,
              WORK.asMatrix(), LDA, RWORK, RCONDC, RESULT(2));

          // Print information about the tests that did not pass
          // the threshold.

          for (var K = 1; K <= 2; K++) {
            if (RESULT[K] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.println(
                  ' UPLO = \'${UPLO.a1}\', N =${N.i5}, NB =${NB.i4}, type ${IMAT.i2}, test ${K.i2}, ratio =${RESULT[K].g12_5}');
              NFAIL++;
            }
          }
          NRUN += 2;

          // Skip the rest of the tests unless this is the first
          // blocksize.

          if (INB != 1) continue;

          for (var IRHS = 1; IRHS <= NNS; IRHS++) {
            final NRHS = NSVAL[IRHS];

            // +    TEST 3
            // Solve and compute residual for A * X = B .

            srnamc.SRNAMT = 'DLARHS';
            dlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A.asMatrix(),
                LDA, XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

            srnamc.SRNAMT = 'DPOTRS';
            dpotrs(
                UPLO, N, NRHS, AFAC.asMatrix(), LDA, X.asMatrix(), LDA, INFO);

            // Check error code from DPOTRS.

            if (INFO.value != 0) {
              alaerh(PATH, 'DPOTRS', INFO.value, 0, UPLO, N, N, -1, -1, NRHS,
                  IMAT, NFAIL, NERRS, NOUT);
            }

            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
            dpot02(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                WORK.asMatrix(), LDA, RWORK, RESULT(3));

            // +    TEST 4
            // Check solution from generated exact solution.

            dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                RCONDC.value, RESULT(4));

            // +    TESTS 5, 6, and 7
            // Use iterative refinement to improve the solution.

            srnamc.SRNAMT = 'DPORFS';
            dporfs(
                UPLO,
                N,
                NRHS,
                A.asMatrix(),
                LDA,
                AFAC.asMatrix(),
                LDA,
                B.asMatrix(),
                LDA,
                X.asMatrix(),
                LDA,
                RWORK,
                RWORK(NRHS + 1),
                WORK,
                IWORK,
                INFO);

            // Check error code from DPORFS.

            if (INFO.value != 0) {
              alaerh(PATH, 'DPORFS', INFO.value, 0, UPLO, N, N, -1, -1, NRHS,
                  IMAT, NFAIL, NERRS, NOUT);
            }

            dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                RCONDC.value, RESULT(5));
            dpot05(
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
                RESULT(6));

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

          // +    TEST 8
          // Get an estimate of RCOND = 1/COND.

          final ANORM = dlansy('1', UPLO, N, A.asMatrix(), LDA, RWORK);
          final RCOND = Box(0.0);
          srnamc.SRNAMT = 'DPOCON';
          dpocon(
              UPLO, N, AFAC.asMatrix(), LDA, ANORM, RCOND, WORK, IWORK, INFO);

          // Check error code from DPOCON.

          if (INFO.value != 0) {
            alaerh(PATH, 'DPOCON', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
                NFAIL, NERRS, NOUT);
          }

          RESULT[8] = dget06(RCOND.value, RCONDC.value);

          // Print the test ratio if it is >= THRESH.

          if (RESULT[8] >= THRESH) {
            if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
            NOUT.println(
                ' UPLO = \'${UPLO.a1}\', N =${N.i5},${' ' * 10} type ${IMAT.i2}, test(${8.i2}) =${RESULT[8].g12_5}');
            NFAIL++;
          }
          NRUN++;
        }
      }
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
