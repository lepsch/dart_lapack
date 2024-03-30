import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlantp.dart';
import 'package:lapack/src/dlatps.dart';
import 'package:lapack/src/dtpcon.dart';
import 'package:lapack/src/dtprfs.dart';
import 'package:lapack/src/dtptri.dart';
import 'package:lapack/src/dtptrs.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'derrtr.dart';
import 'dget04.dart';
import 'dlarhs.dart';
import 'dlattp.dart';
import 'dtpt01.dart';
import 'dtpt02.dart';
import 'dtpt03.dart';
import 'dtpt05.dart';
import 'dtpt06.dart';

void dchktp(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NNS,
  final Array<int> NSVAL_,
  final double THRESH,
  final bool TSTERR,
  final int NMAX,
  final Array<double> AP_,
  final Array<double> AINVP_,
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
  final NSVAL = NSVAL_.having();
  final AP = AP_.having();
  final AINVP = AINVP_.having();
  final B = B_.having();
  final X = X_.having();
  final XACT = XACT_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  const NTYPE1 = 10, NTYPES = 18;
  const NTESTS = 9;
  const NTRAN = 3;
  const ONE = 1.0, ZERO = 0.0;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'], TRANSS = ['N', 'T', 'C'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}TP';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) derrtr(PATH, NOUT);
  infoc.INFOT = 0;

  for (var IN = 1; IN <= NN; IN++) {
    // Do for each value of N in NVAL

    final N = NVAL[IN];
    final LDA = max(1, N);
    final LAP = LDA * (LDA + 1) ~/ 2;
    var XTYPE = 'N';

    for (var IMAT = 1; IMAT <= NTYPE1; IMAT++) {
      // Do the tests only if DOTYPE( IMAT ) is true.

      if (!DOTYPE[IMAT]) continue;

      for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
        // Do first for UPLO = 'U', then for UPLO = 'L'

        final UPLO = UPLOS[IUPLO];

        // Call DLATTP to generate a triangular test matrix.

        final DIAG = Box('');
        srnamc.SRNAMT = 'DLATTP';
        dlattp(IMAT, UPLO, 'No transpose', DIAG, ISEED, N, AP, X, WORK, INFO);

        // Set IDIAG = 1 for non-unit matrices, 2 for unit.

        final IDIAG = lsame(DIAG.value, 'N') ? 1 : 2;

        // +    TEST 1
        // Form the inverse of A.

        if (N > 0) dcopy(LAP, AP, 1, AINVP, 1);
        srnamc.SRNAMT = 'DTPTRI';
        dtptri(UPLO, DIAG.value, N, AINVP, INFO);

        // Check error code from DTPTRI.

        if (INFO.value != 0) {
          alaerh(PATH, 'DTPTRI', INFO.value, 0, UPLO + DIAG.value, N, N, -1, -1,
              -1, IMAT, NFAIL, NERRS, NOUT);
        }

        // Compute the infinity-norm condition number of A.

        final ANORM = dlantp('I', UPLO, DIAG.value, N, AP, RWORK);
        final AINVNM = dlantp('I', UPLO, DIAG.value, N, AINVP, RWORK);
        final double RCONDI;
        if (ANORM <= ZERO || AINVNM <= ZERO) {
          RCONDI = ONE;
        } else {
          RCONDI = (ONE / ANORM) / AINVNM;
        }

        // Compute the residual for the triangular matrix times its
        // inverse.  Also compute the 1-norm condition number of A.

        final RCONDO = Box(ZERO);
        dtpt01(UPLO, DIAG.value, N, AP, AINVP, RCONDO, RWORK, RESULT(1));

        // Print the test ratio if it is >= THRESH.

        if (RESULT[1] >= THRESH) {
          if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
          NOUT.println(
              ' UPLO=\'${UPLO.a1}\', DIAG=\'${DIAG.value.a1}\', N=${N.i5}, type ${IMAT.i2}, test(${1.i2})= ${RESULT[1].g12_5}');
          NFAIL++;
        }
        NRUN++;

        for (var IRHS = 1; IRHS <= NNS; IRHS++) {
          final NRHS = NSVAL[IRHS];
          XTYPE = 'N';

          for (var ITRAN = 1; ITRAN <= NTRAN; ITRAN++) {
            // Do for op(A) = A, A**T, or A**H.

            final TRANS = TRANSS[ITRAN - 1];
            final (_, RCONDC) =
                ITRAN == 1 ? ('O', RCONDO.value) : ('I', RCONDI);

            // +    TEST 2
            // Solve and compute residual for op(A)*x = b.

            srnamc.SRNAMT = 'DLARHS';
            dlarhs(
                PATH,
                XTYPE,
                UPLO,
                TRANS,
                N,
                N,
                0,
                IDIAG,
                NRHS,
                AP.asMatrix(),
                LAP,
                XACT.asMatrix(),
                LDA,
                B.asMatrix(),
                LDA,
                ISEED,
                INFO);
            XTYPE = 'C';
            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

            srnamc.SRNAMT = 'DTPTRS';
            dtptrs(
                UPLO, TRANS, DIAG.value, N, NRHS, AP, X.asMatrix(), LDA, INFO);

            // Check error code from DTPTRS.

            if (INFO.value != 0) {
              alaerh(PATH, 'DTPTRS', INFO.value, 0, UPLO + TRANS + DIAG.value,
                  N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT);
            }

            dtpt02(UPLO, TRANS, DIAG.value, N, NRHS, AP, X.asMatrix(), LDA,
                B.asMatrix(), LDA, WORK, RESULT(2));

            // +    TEST 3
            // Check solution from generated exact solution.

            dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                RESULT(3));

            // +    TESTS 4, 5, and 6
            // Use iterative refinement to improve the solution and
            // compute error bounds.

            srnamc.SRNAMT = 'DTPRFS';
            dtprfs(UPLO, TRANS, DIAG.value, N, NRHS, AP, B.asMatrix(), LDA,
                X.asMatrix(), LDA, RWORK, RWORK(NRHS + 1), WORK, IWORK, INFO);

            // Check error code from DTPRFS.

            if (INFO.value != 0) {
              alaerh(PATH, 'DTPRFS', INFO.value, 0, UPLO + TRANS + DIAG.value,
                  N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT);
            }

            dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                RESULT(4));
            dtpt05(
                UPLO,
                TRANS,
                DIAG.value,
                N,
                NRHS,
                AP,
                B.asMatrix(),
                LDA,
                X.asMatrix(),
                LDA,
                XACT.asMatrix(),
                LDA,
                RWORK,
                RWORK(NRHS + 1),
                RESULT(5));

            // Print information about the tests that did not pass
            // the threshold.

            for (var K = 2; K <= 6; K++) {
              if (RESULT[K] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                NOUT.println(
                    ' UPLO=\'${UPLO.a1}\', TRANS=\'${TRANS.a1}\', DIAG=\'${DIAG.value.a1}\', N=${N.i5}\', NRHS=${NRHS.i5}, type ${IMAT.i2}, test(${K.i2})= ${RESULT[K].g12_5}');
                NFAIL++;
              }
            }
            NRUN += 5;
          }
        }

        // +    TEST 7
        // Get an estimate of RCOND = 1/CNDNUM.

        for (var ITRAN = 1; ITRAN <= 2; ITRAN++) {
          final (NORM, RCONDC) =
              ITRAN == 1 ? ('O', RCONDO.value) : ('I', RCONDI);

          final RCOND = Box(ZERO);
          srnamc.SRNAMT = 'DTPCON';
          dtpcon(NORM, UPLO, DIAG.value, N, AP, RCOND, WORK, IWORK, INFO);

          // Check error code from DTPCON.

          if (INFO.value != 0) {
            alaerh(PATH, 'DTPCON', INFO.value, 0, NORM + UPLO + DIAG.value, N,
                N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT);
          }

          dtpt06(
              RCOND.value, RCONDC, UPLO, DIAG.value, N, AP, RWORK, RESULT(7));

          // Print the test ratio if it is >= THRESH.

          if (RESULT[7] >= THRESH) {
            if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
            NOUT.println(
                ' DTPCON( \'${NORM.a1}\'${UPLO.a1}\'${DIAG.value.a1}\',${N.i5}, ... ), type ${IMAT.i2}, test(${7.i2})=${RESULT[7].g12_5}');
            NFAIL++;
          }
          NRUN++;
        }
      }
    }

    // Use pathological test matrices to test DLATPS.

    for (var IMAT = NTYPE1 + 1; IMAT <= NTYPES; IMAT++) {
      // Do the tests only if DOTYPE( IMAT ) is true.

      if (!DOTYPE[IMAT]) continue;

      for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
        // Do first for UPLO = 'U', then for UPLO = 'L'

        final UPLO = UPLOS[IUPLO - 1];
        for (var ITRAN = 1; ITRAN <= NTRAN; ITRAN++) {
          // Do for op(A) = A, A**T, or A**H.

          final TRANS = TRANSS[ITRAN - 1];

          // Call DLATTP to generate a triangular test matrix.

          final DIAG = Box('');
          srnamc.SRNAMT = 'DLATTP';
          dlattp(IMAT, UPLO, TRANS, DIAG, ISEED, N, AP, X, WORK, INFO);

          // +    TEST 8
          // Solve the system op(A)*x = b.

          srnamc.SRNAMT = 'DLATPS';
          dcopy(N, X, 1, B, 1);
          final SCALE = Box(ZERO);
          dlatps(UPLO, TRANS, DIAG.value, 'N', N, AP, B, SCALE, RWORK, INFO);

          // Check error code from DLATPS.

          if (INFO.value != 0) {
            alaerh(PATH, 'DLATPS', INFO.value, 0, '$UPLO$TRANS${DIAG}N', N, N,
                -1, -1, -1, IMAT, NFAIL, NERRS, NOUT);
          }

          dtpt03(UPLO, TRANS, DIAG.value, N, 1, AP, SCALE.value, RWORK, ONE,
              B.asMatrix(), LDA, X.asMatrix(), LDA, WORK, RESULT(8));

          // +    TEST 9
          // Solve op(A)*x = b again with NORMIN = 'Y'.

          dcopy(N, X, 1, B(N + 1), 1);
          dlatps(UPLO, TRANS, DIAG.value, 'Y', N, AP, B(N + 1), SCALE, RWORK,
              INFO);

          // Check error code from DLATPS.

          if (INFO.value != 0) {
            alaerh(PATH, 'DLATPS', INFO.value, 0, '$UPLO$TRANS${DIAG}Y', N, N,
                -1, -1, -1, IMAT, NFAIL, NERRS, NOUT);
          }

          dtpt03(UPLO, TRANS, DIAG.value, N, 1, AP, SCALE.value, RWORK, ONE,
              B(N + 1).asMatrix(), LDA, X.asMatrix(), LDA, WORK, RESULT(9));

          // Print information about the tests that did not pass
          // the threshold.

          void printFailedTest(String s, String b, int test, double ratio) {
            NOUT.println(
                ' $s( \'${UPLO.a1}\'${TRANS.a1}\'${DIAG.value.a1}\'${b.a1}\',${N.i5}, ... ), type ${IMAT.i2}, test(${test.i2})=${ratio.g12_5}');
          }

          if (RESULT[8] >= THRESH) {
            if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
            printFailedTest('DLATPS', 'N', 8, RESULT[8]);
            NFAIL++;
          }
          if (RESULT[9] >= THRESH) {
            if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
            printFailedTest('DLATPS', 'Y', 9, RESULT[9]);
            NFAIL++;
          }
          NRUN += 2;
        }
      }
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
