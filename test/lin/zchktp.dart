import 'dart:math';

import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlantp.dart';
import 'package:lapack/src/zlatps.dart';
import 'package:lapack/src/ztpcon.dart';
import 'package:lapack/src/ztprfs.dart';
import 'package:lapack/src/ztptri.dart';
import 'package:lapack/src/ztptrs.dart';

import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'zerrtr.dart';
import 'zget04.dart';
import 'zlarhs.dart';
import 'zlattp.dart';
import 'ztpt01.dart';
import 'ztpt02.dart';
import 'ztpt03.dart';
import 'ztpt05.dart';
import 'ztpt06.dart';

void zchktp(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NNS,
  final Array<int> NSVAL_,
  final double THRESH,
  final bool TSTERR,
  final int NMAX,
  final Array<Complex> AP_,
  final Array<Complex> AINVP_,
  final Array<Complex> B_,
  final Array<Complex> X_,
  final Array<Complex> XACT_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
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

  final PATH = '${'Zomplex precision'[0]}TP';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) zerrtr(PATH, NOUT);
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

        final UPLO = UPLOS[IUPLO - 1];

        // Call ZLATTP to generate a triangular test matrix.

        srnamc.SRNAMT = 'ZLATTP';
        final DIAG = Box('');
        zlattp(IMAT, UPLO, 'No transpose', DIAG, ISEED, N, AP, X, WORK, RWORK,
            INFO);

        // Set IDIAG = 1 for non-unit matrices, 2 for unit.

        final IDIAG = lsame(DIAG.value, 'N') ? 1 : 2;

        // +    TEST 1
        // Form the inverse of A.

        if (N > 0) zcopy(LAP, AP, 1, AINVP, 1);
        srnamc.SRNAMT = 'ZTPTRI';
        ztptri(UPLO, DIAG.value, N, AINVP, INFO);

        // Check error code from ZTPTRI.

        if (INFO.value != 0) {
          alaerh(PATH, 'ZTPTRI', INFO.value, 0, UPLO + DIAG.value, N, N, -1, -1,
              -1, IMAT, NFAIL, NERRS, NOUT);
        }

        // Compute the infinity-norm condition number of A.

        final ANORM = zlantp('I', UPLO, DIAG.value, N, AP, RWORK);
        final AINVNM = zlantp('I', UPLO, DIAG.value, N, AINVP, RWORK);
        final double RCONDI;
        if (ANORM <= ZERO || AINVNM <= ZERO) {
          RCONDI = ONE;
        } else {
          RCONDI = (ONE / ANORM) / AINVNM;
        }

        // Compute the residual for the triangular matrix times its
        // inverse.  Also compute the 1-norm condition number of A.

        final RCONDO = Box(ZERO);
        ztpt01(UPLO, DIAG.value, N, AP, AINVP, RCONDO, RWORK, RESULT(1));

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

            srnamc.SRNAMT = 'ZLARHS';
            zlarhs(
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
            zlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

            srnamc.SRNAMT = 'ZTPTRS';
            ztptrs(
                UPLO, TRANS, DIAG.value, N, NRHS, AP, X.asMatrix(), LDA, INFO);

            // Check error code from ZTPTRS.

            if (INFO.value != 0) {
              alaerh(PATH, 'ZTPTRS', INFO.value, 0, UPLO + TRANS + DIAG.value,
                  N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT);
            }

            ztpt02(UPLO, TRANS, DIAG.value, N, NRHS, AP, X.asMatrix(), LDA,
                B.asMatrix(), LDA, WORK, RWORK, RESULT(2));

            // +    TEST 3
            // Check solution from generated exact solution.

            zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                RESULT(3));

            // +    TESTS 4, 5, and 6
            // Use iterative refinement to improve the solution and
            // compute error bounds.

            srnamc.SRNAMT = 'ZTPRFS';
            ztprfs(
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
                RWORK,
                RWORK(NRHS + 1),
                WORK,
                RWORK(2 * NRHS + 1),
                INFO);

            // Check error code from ZTPRFS.

            if (INFO.value != 0) {
              alaerh(PATH, 'ZTPRFS', INFO.value, 0, UPLO + TRANS + DIAG.value,
                  N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT);
            }

            zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                RESULT(4));
            ztpt05(
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

          srnamc.SRNAMT = 'ZTPCON';
          final RCOND = Box(ZERO);
          ztpcon(NORM, UPLO, DIAG.value, N, AP, RCOND, WORK, RWORK, INFO);

          // Check error code from ZTPCON.

          if (INFO.value != 0) {
            alaerh(PATH, 'ZTPCON', INFO.value, 0, NORM + UPLO + DIAG.value, N,
                N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT);
          }

          ztpt06(
              RCOND.value, RCONDC, UPLO, DIAG.value, N, AP, RWORK, RESULT(7));

          // Print the test ratio if it is >= THRESH.

          if (RESULT[7] >= THRESH) {
            if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
            NOUT.println(
                ' ZTPCON( \'${NORM.a1}\'${UPLO.a1}\'${DIAG.value.a1}\',${N.i5}, ... ), type ${IMAT.i2}, test(${7.i2})=${RESULT[7].g12_5}');
            NFAIL++;
          }
          NRUN++;
        }
      }
    }

    // Use pathological test matrices to test ZLATPS.

    for (var IMAT = NTYPE1 + 1; IMAT <= NTYPES; IMAT++) {
      // Do the tests only if DOTYPE( IMAT ) is true.

      if (!DOTYPE[IMAT]) continue;

      for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
        // Do first for UPLO = 'U', then for UPLO = 'L'

        final UPLO = UPLOS[IUPLO - 1];
        for (var ITRAN = 1; ITRAN <= NTRAN; ITRAN++) {
          // Do for op(A) = A, A**T, or A**H.

          final TRANS = TRANSS[ITRAN - 1];

          // Call ZLATTP to generate a triangular test matrix.

          srnamc.SRNAMT = 'ZLATTP';
          final DIAG = Box('');
          zlattp(IMAT, UPLO, TRANS, DIAG, ISEED, N, AP, X, WORK, RWORK, INFO);

          // +    TEST 8
          // Solve the system op(A)*x = b.

          srnamc.SRNAMT = 'ZLATPS';
          zcopy(N, X, 1, B, 1);
          final SCALE = Box(ZERO);
          zlatps(UPLO, TRANS, DIAG.value, 'N', N, AP, B, SCALE, RWORK, INFO);

          // Check error code from ZLATPS.

          if (INFO.value != 0) {
            alaerh(PATH, 'ZLATPS', INFO.value, 0, '$UPLO$TRANS${DIAG.value}N',
                N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT);
          }

          ztpt03(UPLO, TRANS, DIAG.value, N, 1, AP, SCALE.value, RWORK, ONE,
              B.asMatrix(), LDA, X.asMatrix(), LDA, WORK, RESULT(8));

          // +    TEST 9
          // Solve op(A)*x = b again with NORMIN = 'Y'.

          zcopy(N, X, 1, B(N + 1), 1);
          zlatps(UPLO, TRANS, DIAG.value, 'Y', N, AP, B(N + 1), SCALE, RWORK,
              INFO);

          // Check error code from ZLATPS.

          if (INFO.value != 0) {
            alaerh(PATH, 'ZLATPS', INFO.value, 0, '$UPLO$TRANS${DIAG.value}Y',
                N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT);
          }

          ztpt03(UPLO, TRANS, DIAG.value, N, 1, AP, SCALE.value, RWORK, ONE,
              B(N + 1).asMatrix(), LDA, X.asMatrix(), LDA, WORK, RESULT(9));

          // Print information about the tests that did not pass
          // the threshold.

          if (RESULT[8] >= THRESH) {
            if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
            NOUT.println(
                ' ZLATPS( \'${UPLO.a1}\'${TRANS.a1}\'${DIAG.value.a1}\'${'N'.a1}\',${N.i5}, ... ), type ${IMAT.i2}, test(${8.i2})=${RESULT[8].g12_5}');
            NFAIL++;
          }
          if (RESULT[9] >= THRESH) {
            if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
            NOUT.println(
                ' ZLATPS( \'${UPLO.a1}\'${TRANS.a1}\'${DIAG.value.a1}\'${'Y'.a1}\',${N.i5}, ... ), type ${IMAT.i2}, test(${9.i2})=${RESULT[9].g12_5}');
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
