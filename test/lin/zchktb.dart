// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zcopy.dart';
import 'package:dart_lapack/src/blas/ztbsv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlantb.dart';
import 'package:dart_lapack/src/zlantr.dart';
import 'package:dart_lapack/src/zlaset.dart';
import 'package:dart_lapack/src/zlatbs.dart';
import 'package:dart_lapack/src/ztbcon.dart';
import 'package:dart_lapack/src/ztbrfs.dart';
import 'package:dart_lapack/src/ztbtrs.dart';

import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'zerrtr.dart';
import 'zget04.dart';
import 'zlarhs.dart';
import 'zlattb.dart';
import 'ztbt02.dart';
import 'ztbt03.dart';
import 'ztbt05.dart';
import 'ztbt06.dart';

void zchktb(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NNS,
  final Array<int> NSVAL_,
  final double THRESH,
  final bool TSTERR,
  final int NMAX,
  final Array<Complex> AB_,
  final Array<Complex> AINV_,
  final Array<Complex> B_,
  final Array<Complex> X_,
  final Array<Complex> XACT_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Nout NOUT,
) {
  final DOTYPE = DOTYPE_.having();
  final NVAL = NVAL_.having();
  final NSVAL = NSVAL_.having();
  final AB = AB_.having();
  final AINV = AINV_.having();
  final B = B_.having();
  final X = X_.having();
  final XACT = XACT_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const NTYPE1 = 9, NTYPES = 17, NTESTS = 8, NTRAN = 3;
  const ONE = 1.0, ZERO = 0.0;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'], TRANSS = ['N', 'T', 'C'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Zomplex precision'[0]}TB';
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
    var XTYPE = 'N';
    final NIMAT = N <= 0 ? 1 : NTYPE1;
    final NIMAT2 = N <= 0 ? NTYPE1 + 1 : NTYPES;

    final NK = min(N + 1, 4);
    for (var IK = 1; IK <= NK; IK++) {
      // Do for KD = 0, N, (3N-1)/4, and (N+1)/4. This order makes
      // it easier to skip redundant values for small values of N.

      final KD = switch (IK) {
        1 => 0,
        2 => max(N, 0),
        3 => (3 * N - 1) ~/ 4,
        4 => (N + 1) ~/ 4,
        _ => throw UnimplementedError(),
      };
      final LDAB = KD + 1;

      for (var IMAT = 1; IMAT <= NIMAT; IMAT++) {
        // Do the tests only if DOTYPE( IMAT ) is true.

        if (!DOTYPE[IMAT]) continue;

        for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
          // Do first for UPLO = 'U', then for UPLO = 'L'

          final UPLO = UPLOS[IUPLO - 1];

          // Call ZLATTB to generate a triangular test matrix.

          srnamc.SRNAMT = 'ZLATTB';
          final DIAG = Box('');
          zlattb(IMAT, UPLO, 'No transpose', DIAG, ISEED, N, KD, AB.asMatrix(),
              LDAB, X, WORK, RWORK, INFO);

          // Set IDIAG = 1 for non-unit matrices, 2 for unit.

          final IDIAG = lsame(DIAG.value, 'N') ? 1 : 2;

          // Form the inverse of A so we can get a good estimate
          // of RCONDC = 1/(norm(A) * norm(inv(A))).

          zlaset('Full', N, N, Complex.zero, Complex.one, AINV.asMatrix(), LDA);
          if (lsame(UPLO, 'U')) {
            for (var J = 1; J <= N; J++) {
              ztbsv(UPLO, 'No transpose', DIAG.value, J, KD, AB.asMatrix(),
                  LDAB, AINV((J - 1) * LDA + 1), 1);
            }
          } else {
            for (var J = 1; J <= N; J++) {
              ztbsv(
                  UPLO,
                  'No transpose',
                  DIAG.value,
                  N - J + 1,
                  KD,
                  AB((J - 1) * LDAB + 1).asMatrix(),
                  LDAB,
                  AINV((J - 1) * LDA + J),
                  1);
            }
          }

          // Compute the 1-norm condition number of A.

          var ANORM =
              zlantb('1', UPLO, DIAG.value, N, KD, AB.asMatrix(), LDAB, RWORK);
          var AINVNM =
              zlantr('1', UPLO, DIAG.value, N, N, AINV.asMatrix(), LDA, RWORK);
          final double RCONDO;
          if (ANORM <= ZERO || AINVNM <= ZERO) {
            RCONDO = ONE;
          } else {
            RCONDO = (ONE / ANORM) / AINVNM;
          }

          // Compute the infinity-norm condition number of A.

          ANORM =
              zlantb('I', UPLO, DIAG.value, N, KD, AB.asMatrix(), LDAB, RWORK);
          AINVNM =
              zlantr('I', UPLO, DIAG.value, N, N, AINV.asMatrix(), LDA, RWORK);
          final double RCONDI;
          if (ANORM <= ZERO || AINVNM <= ZERO) {
            RCONDI = ONE;
          } else {
            RCONDI = (ONE / ANORM) / AINVNM;
          }

          for (var IRHS = 1; IRHS <= NNS; IRHS++) {
            final NRHS = NSVAL[IRHS];
            XTYPE = 'N';

            for (var ITRAN = 1; ITRAN <= NTRAN; ITRAN++) {
              // Do for op(A) = A, A**T, or A**H.

              final TRANS = TRANSS[ITRAN - 1];
              final (_, RCONDC) = ITRAN == 1 ? ('O', RCONDO) : ('I', RCONDI);

              // +    TEST 1
              // Solve and compute residual for op(A)*x = b.

              srnamc.SRNAMT = 'ZLARHS';
              zlarhs(
                  PATH,
                  XTYPE,
                  UPLO,
                  TRANS,
                  N,
                  N,
                  KD,
                  IDIAG,
                  NRHS,
                  AB.asMatrix(),
                  LDAB,
                  XACT.asMatrix(),
                  LDA,
                  B.asMatrix(),
                  LDA,
                  ISEED,
                  INFO);
              XTYPE = 'C';
              zlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

              srnamc.SRNAMT = 'ZTBTRS';
              ztbtrs(UPLO, TRANS, DIAG.value, N, KD, NRHS, AB.asMatrix(), LDAB,
                  X.asMatrix(), LDA, INFO);

              // Check error code from ZTBTRS.

              if (INFO.value != 0) {
                alaerh(PATH, 'ZTBTRS', INFO.value, 0, UPLO + TRANS + DIAG.value,
                    N, N, KD, KD, NRHS, IMAT, NFAIL, NERRS, NOUT);
              }

              ztbt02(UPLO, TRANS, DIAG.value, N, KD, NRHS, AB.asMatrix(), LDAB,
                  X.asMatrix(), LDA, B.asMatrix(), LDA, WORK, RWORK, RESULT(1));

              // +    TEST 2
              // Check solution from generated exact solution.

              zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                  RESULT(2));

              // +    TESTS 3, 4, and 5
              // Use iterative refinement to improve the solution
              // and compute error bounds.

              srnamc.SRNAMT = 'ZTBRFS';
              ztbrfs(
                  UPLO,
                  TRANS,
                  DIAG.value,
                  N,
                  KD,
                  NRHS,
                  AB.asMatrix(),
                  LDAB,
                  B.asMatrix(),
                  LDA,
                  X.asMatrix(),
                  LDA,
                  RWORK,
                  RWORK(NRHS + 1),
                  WORK,
                  RWORK(2 * NRHS + 1),
                  INFO);

              // Check error code from ZTBRFS.

              if (INFO.value != 0) {
                alaerh(PATH, 'ZTBRFS', INFO.value, 0, UPLO + TRANS + DIAG.value,
                    N, N, KD, KD, NRHS, IMAT, NFAIL, NERRS, NOUT);
              }

              zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                  RESULT(3));
              ztbt05(
                  UPLO,
                  TRANS,
                  DIAG.value,
                  N,
                  KD,
                  NRHS,
                  AB.asMatrix(),
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

              // Print information about the tests that did not
              // pass the threshold.

              for (var K = 1; K <= 5; K++) {
                if (RESULT[K] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                  NOUT.println(
                      ' UPLO=\'${UPLO.a1}\', TRANS=\'${TRANS.a1}\', DIAG=\'${DIAG.value.a1}\', N=${N.i5}, KD=${KD.i5}, NRHS=${NRHS.i5}, type ${IMAT.i2}, test(${K.i2})=${RESULT[K].g12_5}');
                  NFAIL++;
                }
              }
              NRUN += 5;
            }
          }

          // +    TEST 6
          // Get an estimate of RCOND = 1/CNDNUM.

          for (var ITRAN = 1; ITRAN <= 2; ITRAN++) {
            final (NORM, RCONDC) = ITRAN == 1 ? ('O', RCONDO) : ('I', RCONDI);

            srnamc.SRNAMT = 'ZTBCON';
            final RCOND = Box(ZERO);
            ztbcon(NORM, UPLO, DIAG.value, N, KD, AB.asMatrix(), LDAB, RCOND,
                WORK, RWORK, INFO);

            // Check error code from ZTBCON.

            if (INFO.value != 0) {
              alaerh(PATH, 'ZTBCON', INFO.value, 0, NORM + UPLO + DIAG.value, N,
                  N, KD, KD, -1, IMAT, NFAIL, NERRS, NOUT);
            }

            ztbt06(RCOND.value, RCONDC, UPLO, DIAG.value, N, KD, AB.asMatrix(),
                LDAB, RWORK, RESULT(6));

            // Print the test ratio if it is >= THRESH.

            if (RESULT[6] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.println(
                  ' ZTBCON( \'${NORM.a1}\'${UPLO.a1}\'${DIAG.value.a1}\',${N.i5},${KD.i5},  ... ), type ${IMAT.i2}, test(${6.i2})=${RESULT[6].g12_5}');
              NFAIL++;
            }
            NRUN++;
          }
        }
      }

      // Use pathological test matrices to test ZLATBS.

      for (var IMAT = NTYPE1 + 1; IMAT <= NIMAT2; IMAT++) {
        // Do the tests only if DOTYPE( IMAT ) is true.

        if (!DOTYPE[IMAT]) continue;

        for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
          // Do first for UPLO = 'U', then for UPLO = 'L'

          final UPLO = UPLOS[IUPLO - 1];
          for (var ITRAN = 1; ITRAN <= NTRAN; ITRAN++) {
            // Do for op(A) = A, A**T, and A**H.

            final TRANS = TRANSS[ITRAN - 1];

            // Call ZLATTB to generate a triangular test matrix.

            srnamc.SRNAMT = 'ZLATTB';
            final DIAG = Box('');
            zlattb(IMAT, UPLO, TRANS, DIAG, ISEED, N, KD, AB.asMatrix(), LDAB,
                X, WORK, RWORK, INFO);

            // +    TEST 7
            // Solve the system op(A)*x = b

            srnamc.SRNAMT = 'ZLATBS';
            zcopy(N, X, 1, B, 1);
            final SCALE = Box(ZERO);
            zlatbs(UPLO, TRANS, DIAG.value, 'N', N, KD, AB.asMatrix(), LDAB, B,
                SCALE, RWORK, INFO);

            // Check error code from ZLATBS.

            if (INFO.value != 0) {
              alaerh(PATH, 'ZLATBS', INFO.value, 0, '$UPLO$TRANS${DIAG.value}N',
                  N, N, KD, KD, -1, IMAT, NFAIL, NERRS, NOUT);
            }

            ztbt03(
                UPLO,
                TRANS,
                DIAG.value,
                N,
                KD,
                1,
                AB.asMatrix(),
                LDAB,
                SCALE.value,
                RWORK,
                ONE,
                B.asMatrix(),
                LDA,
                X.asMatrix(),
                LDA,
                WORK,
                RESULT(7));

            // +    TEST 8
            // Solve op(A)*x = b again with NORMIN = 'Y'.

            zcopy(N, X, 1, B, 1);
            zlatbs(UPLO, TRANS, DIAG.value, 'Y', N, KD, AB.asMatrix(), LDAB, B,
                SCALE, RWORK, INFO);

            // Check error code from ZLATBS.

            if (INFO.value != 0) {
              alaerh(PATH, 'ZLATBS', INFO.value, 0, '$UPLO$TRANS${DIAG.value}Y',
                  N, N, KD, KD, -1, IMAT, NFAIL, NERRS, NOUT);
            }

            ztbt03(
                UPLO,
                TRANS,
                DIAG.value,
                N,
                KD,
                1,
                AB.asMatrix(),
                LDAB,
                SCALE.value,
                RWORK,
                ONE,
                B.asMatrix(),
                LDA,
                X.asMatrix(),
                LDA,
                WORK,
                RESULT(8));

            // Print information about the tests that did not pass
            // the threshold.

            if (RESULT[7] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.println(
                  ' ZLATBS( \'${UPLO.a1}\'${TRANS.a1}\'${DIAG.value.a1}\'${'N'.a1}\',${N.i5},${KD.i5}, ...  ),  type ${IMAT.i2}, test(${7.i1})=${RESULT[7].g12_5}');
              NFAIL++;
            }
            if (RESULT[8] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.println(
                  ' ZLATBS( \'${UPLO.a1}\'${TRANS.a1}\'${DIAG.value.a1}\'${'Y'.a1}\',${N.i5},${KD.i5}, ...  ),  type ${IMAT.i2}, test(${8.i1})=${RESULT[8].g12_5}');
              NFAIL++;
            }
            NRUN += 2;
          }
        }
      }
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
