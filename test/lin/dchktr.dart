// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';
import 'package:test/test.dart';

import '../test_driver.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'derrtr.dart';
import 'dget04.dart';
import 'dlarhs.dart';
import 'dlattr.dart';
import 'dtrt01.dart';
import 'dtrt02.dart';
import 'dtrt03.dart';
import 'dtrt05.dart';
import 'dtrt06.dart';
import 'xlaenv.dart';

void dchktr(
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
  final Array<double> AINV_,
  final Array<double> B_,
  final Array<double> X_,
  final Array<double> XACT_,
  final Array<double> WORK_,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Nout NOUT,
  final TestDriver test,
) {
  final DOTYPE = DOTYPE_.having();
  final NVAL = NVAL_.having();
  final NBVAL = NBVAL_.having();
  final NSVAL = NSVAL_.having();
  final A = A_.having();
  final AINV = AINV_.having();
  final B = B_.having();
  final X = X_.having();
  final XACT = XACT_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  const NTYPE1 = 10, NTYPES = 18, NTESTS = 10, NTRAN = 3;
  const ONE = 1.0, ZERO = 0.0;
  final RESULT = Array<double>(NTESTS), SCALE3 = Array<double>(2);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'], TRANSS = ['N', 'T', 'C'];

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}TR';
  final BIGNUM = dlamch('Overflow') / dlamch('Precision');
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  final ISEED = Array.fromList(ISEEDY);

  test.group('error exits', () {
    // Test the error exits
    if (TSTERR) derrtr(PATH, NOUT, test);
    test.tearDown(() {
      infoc.INFOT = 0;
    });
  });

  test.setUp(() {
    xlaenv(2, 2);
  });

  // Do for each value of N in NVAL
  for (final IN in 1.through(NN)) {
    final N = NVAL[IN];
    final LDA = max(1, N);

    for (final IMAT in 1.through(NTYPE1)) {
      // Do the tests only if DOTYPE( IMAT ) is true.
      final skip = !DOTYPE[IMAT];

      test('DCHKTR - 1 (IN=$IN IMAT=$IMAT)', () {
        final INFO = Box(0);

        for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
          // Do first for UPLO = 'U', then for UPLO = 'L'
          final UPLO = UPLOS[IUPLO - 1];

          // Call DLATTR to generate a triangular test matrix.
          final DIAG = Box('');
          srnamc.SRNAMT = 'DLATTR';
          dlattr(IMAT, UPLO, 'No transpose', DIAG, ISEED, N, A.asMatrix(), LDA,
              X, WORK, INFO);

          // Set IDIAG = 1 for non-unit matrices, 2 for unit.
          final IDIAG = lsame(DIAG.value, 'N') ? 1 : 2;

          // Do for each blocksize in NBVAL
          for (var INB = 1; INB <= NNB; INB++) {
            final NB = NBVAL[INB];
            xlaenv(1, NB);

            // +    TEST 1

            // Form the inverse of A.
            dlacpy(UPLO, N, N, A.asMatrix(), LDA, AINV.asMatrix(), LDA);
            srnamc.SRNAMT = 'DTRTRI';
            dtrtri(UPLO, DIAG.value, N, AINV.asMatrix(), LDA, INFO);

            // Check error code from DTRTRI.
            test.expect(INFO.value, 0);
            if (INFO.value != 0) {
              alaerh(PATH, 'DTRTRI', INFO.value, 0, UPLO + DIAG.value, N, N, -1,
                  -1, NB, IMAT, NFAIL, NERRS, NOUT);
            }

            // Compute the infinity-norm condition number of A.
            final ANORM =
                dlantr('I', UPLO, DIAG.value, N, N, A.asMatrix(), LDA, RWORK);
            final AINVNM = dlantr(
                'I', UPLO, DIAG.value, N, N, AINV.asMatrix(), LDA, RWORK);
            final double RCONDI;
            if (ANORM <= ZERO || AINVNM <= ZERO) {
              RCONDI = ONE;
            } else {
              RCONDI = (ONE / ANORM) / AINVNM;
            }

            // Compute the residual for the triangular matrix times
            // its inverse.  Also compute the 1-norm condition number
            // of A.
            final RCONDO = Box(0.0);
            dtrt01(UPLO, DIAG.value, N, A.asMatrix(), LDA, AINV.asMatrix(), LDA,
                RCONDO, RWORK, RESULT(1));

            // Print the test ratio if it is >= THRESH.
            final reason =
                ' UPLO=\'${UPLO.a1}\', DIAG=\'${DIAG.value.a1}\', N=${N.i5}, NB=${NB.i4}, type ${IMAT.i2}, test(${1.i2})= ${RESULT[1].g12_5}';
            test.expect(RESULT[1], lessThan(THRESH), reason: reason);
            if (RESULT[1] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.println(reason);
              NFAIL++;
            }
            NRUN++;

            // Skip remaining tests if not the first block size.
            if (INB != 1) continue;

            for (var IRHS = 1; IRHS <= NNS; IRHS++) {
              final NRHS = NSVAL[IRHS];
              String XTYPE = 'N';

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
                    A.asMatrix(),
                    LDA,
                    XACT.asMatrix(),
                    LDA,
                    B.asMatrix(),
                    LDA,
                    ISEED,
                    INFO);
                XTYPE = 'C';
                dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

                srnamc.SRNAMT = 'DTRTRS';
                dtrtrs(UPLO, TRANS, DIAG.value, N, NRHS, A.asMatrix(), LDA,
                    X.asMatrix(), LDA, INFO);

                // Check error code from DTRTRS.
                test.expect(INFO.value, 0);
                if (INFO.value != 0) {
                  alaerh(
                      PATH,
                      'DTRTRS',
                      INFO.value,
                      0,
                      UPLO + TRANS + DIAG.value,
                      N,
                      N,
                      -1,
                      -1,
                      NRHS,
                      IMAT,
                      NFAIL,
                      NERRS,
                      NOUT);
                }

                dtrt02(UPLO, TRANS, DIAG.value, N, NRHS, A.asMatrix(), LDA,
                    X.asMatrix(), LDA, B.asMatrix(), LDA, WORK, RESULT(2));

                // +    TEST 3
                // Check solution from generated exact solution.

                dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                    RESULT(3));

                // +    TESTS 4, 5, and 6
                // Use iterative refinement to improve the solution
                // and compute error bounds.

                srnamc.SRNAMT = 'DTRRFS';
                dtrrfs(
                    UPLO,
                    TRANS,
                    DIAG.value,
                    N,
                    NRHS,
                    A.asMatrix(),
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

                // Check error code from DTRRFS.
                test.expect(INFO.value, 0);
                if (INFO.value != 0) {
                  alaerh(
                      PATH,
                      'DTRRFS',
                      INFO.value,
                      0,
                      UPLO + TRANS + DIAG.value,
                      N,
                      N,
                      -1,
                      -1,
                      NRHS,
                      IMAT,
                      NFAIL,
                      NERRS,
                      NOUT);
                }

                dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                    RESULT(4));
                dtrt05(
                    UPLO,
                    TRANS,
                    DIAG.value,
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
                    RESULT(5));

                // Print information about the tests that did not
                // pass the threshold.
                for (var K = 2; K <= 6; K++) {
                  final reason =
                      ' UPLO=\'${UPLO.a1}\', TRANS=\'${TRANS.a1}\', DIAG=\'${DIAG.value.a1}\', N=${N.i5}, NB=${NB.i4}, type ${IMAT.i2}, test(${K.i2})= ${RESULT[K].g12_5}';
                  test.expect(RESULT[K], lessThan(THRESH), reason: reason);
                  if (RESULT[K] >= THRESH) {
                    if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                    NOUT.println(reason);
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
              srnamc.SRNAMT = 'DTRCON';
              dtrcon(NORM, UPLO, DIAG.value, N, A.asMatrix(), LDA, RCOND, WORK,
                  IWORK, INFO);

              // Check error code from DTRCON.
              test.expect(INFO.value, 0);
              if (INFO.value != 0) {
                alaerh(PATH, 'DTRCON', INFO.value, 0, NORM + UPLO + DIAG.value,
                    N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT);
              }

              dtrt06(RCOND.value, RCONDC, UPLO, DIAG.value, N, A.asMatrix(),
                  LDA, RWORK, RESULT(7));

              // Print the test ratio if it is >= THRESH.
              test.expect(RESULT[7], lessThan(THRESH), reason: reason);
              if (RESULT[7] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                NOUT.println(
                    ' NORM=\'${NORM.a1}\', UPLO =\'${UPLO.a1}\', N=${N.i5},${' ' * 11} type ${IMAT.i2}, test(${7.i2})=${RESULT[7].g12_5}');
                NFAIL++;
              }
              NRUN++;
            }
          }
        }
      }, skip: skip);
    }

    // Use pathological test matrices to test DLATRS.
    for (final IMAT in (NTYPE1 + 1).through(NTYPES)) {
      // Do the tests only if DOTYPE( IMAT ) is true.
      final skip = !DOTYPE[IMAT];

      test('DCHKTR - 2 (IN=$IN IMAT=$IMAT)', () {
        final INFO = Box(0);

        for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
          // Do first for UPLO = 'U', then for UPLO = 'L'
          final UPLO = UPLOS[IUPLO - 1];
          for (var ITRAN = 1; ITRAN <= NTRAN; ITRAN++) {
            // Do for op(A) = A, A**T, and A**H.

            final TRANS = TRANSS[ITRAN - 1];

            // Call DLATTR to generate a triangular test matrix.
            final DIAG = Box('');
            srnamc.SRNAMT = 'DLATTR';
            dlattr(IMAT, UPLO, TRANS, DIAG, ISEED, N, A.asMatrix(), LDA, X,
                WORK, INFO);

            // +    TEST 8
            // Solve the system op(A)*x = b.

            srnamc.SRNAMT = 'DLATRS';
            dcopy(N, X, 1, B, 1);
            final SCALE = Box(0.0);
            dlatrs(UPLO, TRANS, DIAG.value, 'N', N, A.asMatrix(), LDA, B, SCALE,
                RWORK, INFO);

            // Check error code from DLATRS.
            test.expect(INFO.value, 0);
            if (INFO.value != 0) {
              alaerh(PATH, 'DLATRS', INFO.value, 0, '$UPLO$TRANS${DIAG}N', N, N,
                  -1, -1, -1, IMAT, NFAIL, NERRS, NOUT);
            }

            dtrt03(
                UPLO,
                TRANS,
                DIAG.value,
                N,
                1,
                A.asMatrix(),
                LDA,
                SCALE.value,
                RWORK,
                ONE,
                B.asMatrix(),
                LDA,
                X.asMatrix(),
                LDA,
                WORK,
                RESULT(8));

            // +    TEST 9
            // Solve op(A)*X = b again with NORMIN = 'Y'.

            dcopy(N, X, 1, B(N + 1), 1);
            dlatrs(UPLO, TRANS, DIAG.value, 'Y', N, A.asMatrix(), LDA, B(N + 1),
                SCALE, RWORK, INFO);

            // Check error code from DLATRS.
            test.expect(INFO.value, 0);
            if (INFO.value != 0) {
              alaerh(PATH, 'DLATRS', INFO.value, 0, '$UPLO$TRANS${DIAG}Y', N, N,
                  -1, -1, -1, IMAT, NFAIL, NERRS, NOUT);
            }

            dtrt03(
                UPLO,
                TRANS,
                DIAG.value,
                N,
                1,
                A.asMatrix(),
                LDA,
                SCALE.value,
                RWORK,
                ONE,
                B(N + 1).asMatrix(),
                LDA,
                X.asMatrix(),
                LDA,
                WORK,
                RESULT(9));

            // +    TEST 10
            // Solve op(A)*X = B

            srnamc.SRNAMT = 'DLATRS3';
            dcopy(N, X, 1, B, 1);
            dcopy(N, X, 1, B(N + 1), 1);
            dscal(N, BIGNUM, B(N + 1), 1);
            dlatrs3(UPLO, TRANS, DIAG.value, 'N', N, 2, A.asMatrix(), LDA,
                B.asMatrix(), max(1, N), SCALE3, RWORK, WORK, NMAX, INFO);

            // Check error code from DLATRS3.
            test.expect(INFO.value, 0);
            if (INFO.value != 0) {
              alaerh(PATH, 'DLATRS3', INFO.value, 0, '$UPLO$TRANS${DIAG}N', N,
                  N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT);
            }
            dtrt03(
                UPLO,
                TRANS,
                DIAG.value,
                N,
                1,
                A.asMatrix(),
                LDA,
                SCALE3[1],
                RWORK,
                ONE,
                B(1).asMatrix(),
                LDA,
                X.asMatrix(),
                LDA,
                WORK,
                RESULT(10));
            dscal(N, BIGNUM, X, 1);
            final RES = Box(ZERO);
            dtrt03(
                UPLO,
                TRANS,
                DIAG.value,
                N,
                1,
                A.asMatrix(),
                LDA,
                SCALE3[2],
                RWORK,
                ONE,
                B(N + 1).asMatrix(),
                LDA,
                X.asMatrix(),
                LDA,
                WORK,
                RES);
            RESULT[10] = max(RESULT[10], RES.value);

            // Print information about the tests that did not pass
            // the threshold.
            for (final (K, SUBNAM, NORMIN) in [
              (8, 'DLATRS', 'N'),
              (9, 'DLATRS', 'Y'),
              (10, 'DLATRS3', 'N')
            ]) {
              final reason =
                  ' $SUBNAM( \'${UPLO.a1}\'${TRANS.a1}\'${DIAG.value.a1}\'${NORMIN.a1}\',${N.i5}, ... ), type ${IMAT.i2}, test(${K.i2})=${RESULT[K].g12_5}';
              test.expect(RESULT[K], lessThan(THRESH), reason: reason);
              if (RESULT[K] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                NOUT.println(reason);
                NFAIL++;
              }
            }
            NRUN += 3;
          }
        }
      }, skip: skip);
    }
  }

  // Print a summary of the results.
  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
