// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';

import '../matgen/zlatms.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'dget06.dart';
import 'zerrge.dart';
import 'zget04.dart';
import 'zgtt01.dart';
import 'zgtt02.dart';
import 'zgtt05.dart';
import 'zlatb4.dart';

void zchkgt(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NNS,
  final Array<int> NSVAL_,
  final double THRESH,
  final bool TSTERR,
  final Array<Complex> A_,
  final Array<Complex> AF_,
  final Array<Complex> B_,
  final Array<Complex> X_,
  final Array<Complex> XACT_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Nout NOUT,
) {
  final DOTYPE = DOTYPE_.having();
  final NVAL = NVAL_.having();
  final NSVAL = NSVAL_.having();
  final A = A_.having();
  final AF = AF_.having();
  final B = B_.having();
  final X = X_.having();
  final XACT = XACT_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  const NTYPES = 12;
  const NTESTS = 7;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  final Z = Array<Complex>(3);
  const ISEEDY = [0, 0, 0, 1], TRANSS = ['N', 'T', 'C'];
  final INFO = Box(0);

  final PATH = '${'Zomplex precision'[0]}GT';
  var NRUN = 0;
  var NFAIL = 0;
  var NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) zerrge(PATH, NOUT);
  infoc.INFOT = 0;

  for (var IN = 1; IN <= NN; IN++) {
    // Do for each value of N in NVAL.

    final N = NVAL[IN];
    final M = max(N - 1, 0);
    final LDA = max(1, N);
    final NIMAT = N <= 0 ? 1 : NTYPES;

    var IZERO = 0;

    for (var IMAT = 1; IMAT <= NIMAT; IMAT++) {
      // Do the tests only if DOTYPE( IMAT ) is true.

      if (!DOTYPE[IMAT]) continue;

      // Set up parameters with ZLATB4.

      final (:TYPE, :KL, :KU, :ANORM, :MODE, CNDNUM: COND, :DIST) =
          zlatb4(PATH, IMAT, N, N);

      final ZEROT = IMAT >= 8 && IMAT <= 10;
      if (IMAT <= 6) {
        // Types 1-6:  generate matrices of known condition number.

        final KOFF = max(2 - KU, 3 - max(1, N)).toInt();
        srnamc.SRNAMT = 'ZLATMS';
        zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, COND, ANORM, KL, KU, 'Z',
            AF(KOFF).asMatrix(), 3, WORK, INFO);

        // Check the error code from ZLATMS.

        if (INFO.value != 0) {
          alaerh(PATH, 'ZLATMS', INFO.value, 0, ' ', N, N, KL, KU, -1, IMAT,
              NFAIL, NERRS, NOUT);
          continue;
        }
        IZERO = 0;

        if (N > 1) {
          zcopy(N - 1, AF(4), 3, A, 1);
          zcopy(N - 1, AF(3), 3, A(N + M + 1), 1);
        }
        zcopy(N, AF(2), 3, A(M + 1), 1);
      } else {
        // Types 7-12:  generate tridiagonal matrices with
        // unknown condition numbers.

        if (!ZEROT || !DOTYPE[7]) {
          // Generate a matrix with elements whose real and
          // imaginary parts are from [-1,1].

          zlarnv(2, ISEED, N + 2 * M, A);
          if (ANORM != ONE) zdscal(N + 2 * M, ANORM, A, 1);
        } else if (IZERO > 0) {
          // Reuse the last matrix by copying back the zeroed out
          // elements.

          if (IZERO == 1) {
            A[N] = Z[2];
            if (N > 1) A[1] = Z[3];
          } else if (IZERO == N) {
            A[3 * N - 2] = Z[1];
            A[2 * N - 1] = Z[2];
          } else {
            A[2 * N - 2 + IZERO] = Z[1];
            A[N - 1 + IZERO] = Z[2];
            A[IZERO] = Z[3];
          }
        }

        // If IMAT > 7, set one column of the matrix to 0.

        if (!ZEROT) {
          IZERO = 0;
        } else if (IMAT == 8) {
          IZERO = 1;
          Z[2] = A[N];
          A[N] = Complex.zero;
          if (N > 1) {
            Z[3] = A[1];
            A[1] = Complex.zero;
          }
        } else if (IMAT == 9) {
          IZERO = N;
          Z[1] = A[3 * N - 2];
          Z[2] = A[2 * N - 1];
          A[3 * N - 2] = Complex.zero;
          A[2 * N - 1] = Complex.zero;
        } else {
          IZERO = (N + 1) ~/ 2;
          for (var I = IZERO; I <= N - 1; I++) {
            A[2 * N - 2 + I] = Complex.zero;
            A[N - 1 + I] = Complex.zero;
            A[I] = Complex.zero;
          }
          A[3 * N - 2] = Complex.zero;
          A[2 * N - 1] = Complex.zero;
        }
      }

      // +    TEST 1
      // Factor A as L*U and compute the ratio
      //    norm(L*U - A) / (n * norm(A) * EPS )

      zcopy(N + 2 * M, A, 1, AF, 1);
      srnamc.SRNAMT = 'ZGTTRF';
      zgttrf(N, AF, AF(M + 1), AF(N + M + 1), AF(N + 2 * M + 1), IWORK, INFO);

      // Check error code from ZGTTRF.

      if (INFO.value != IZERO) {
        alaerh(PATH, 'ZGTTRF', INFO.value, IZERO, ' ', N, N, 1, 1, -1, IMAT,
            NFAIL, NERRS, NOUT);
      }
      final TRFCON = INFO.value != 0;

      zgtt01(N, A, A(M + 1), A(N + M + 1), AF, AF(M + 1), AF(N + M + 1),
          AF(N + 2 * M + 1), IWORK, WORK.asMatrix(), LDA, RWORK, RESULT(1));

      // Print the test ratio if it is >= THRESH.

      if (RESULT[1] >= THRESH) {
        if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
        NOUT.println(
            '${' ' * 12}N =${N.i5},${' ' * 10} type ${IMAT.i2}, test(${1.i2}) = ${RESULT[1].g12_5}');
        NFAIL++;
      }
      NRUN++;

      double RCONDO = ZERO, RCONDI = ZERO;
      for (var ITRAN = 1; ITRAN <= 2; ITRAN++) {
        final TRANS = TRANSS[ITRAN - 1];
        final NORM = ITRAN == 1 ? 'O' : 'I';
        final ANORM = zlangt(NORM, N, A, A(M + 1), A(N + M + 1));

        final double RCONDC;
        if (!TRFCON) {
          // Use ZGTTRS to solve for one column at a time of
          // inv(A), computing the maximum column sum as we go.

          var AINVNM = ZERO;
          for (var I = 1; I <= N; I++) {
            for (var J = 1; J <= N; J++) {
              X[J] = Complex.zero;
            }
            X[I] = Complex.one;
            zgttrs(TRANS, N, 1, AF, AF(M + 1), AF(N + M + 1), AF(N + 2 * M + 1),
                IWORK, X.asMatrix(), LDA, INFO);
            AINVNM = max(AINVNM, dzasum(N, X, 1));
          }

          // Compute RCONDC = 1 / (norm(A) * norm(inv(A))

          if (ANORM <= ZERO || AINVNM <= ZERO) {
            RCONDC = ONE;
          } else {
            RCONDC = (ONE / ANORM) / AINVNM;
          }
          if (ITRAN == 1) {
            RCONDO = RCONDC;
          } else {
            RCONDI = RCONDC;
          }
        } else {
          RCONDC = ZERO;
        }

        // +    TEST 7
        // Estimate the reciprocal of the condition number of the
        // matrix.

        final RCOND = Box(ZERO);
        srnamc.SRNAMT = 'ZGTCON';
        zgtcon(NORM, N, AF, AF(M + 1), AF(N + M + 1), AF(N + 2 * M + 1), IWORK,
            ANORM, RCOND, WORK, INFO);

        // Check error code from ZGTCON.

        if (INFO.value != 0) {
          alaerh(PATH, 'ZGTCON', INFO.value, 0, NORM, N, N, -1, -1, -1, IMAT,
              NFAIL, NERRS, NOUT);
        }

        RESULT[7] = dget06(RCOND.value, RCONDC);

        // Print the test ratio if it is >= THRESH.

        if (RESULT[7] >= THRESH) {
          if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
          NOUT.println(
              ' NORM =\'${NORM.a1}\', N =${N.i5},${' ' * 10} type ${IMAT.i2}, test(${7.i2}) = ${RESULT[7].g12_5}');
          NFAIL++;
        }
        NRUN++;
      }

      // Skip the remaining tests if the matrix is singular.

      if (TRFCON) continue;

      for (var IRHS = 1; IRHS <= NNS; IRHS++) {
        final NRHS = NSVAL[IRHS];

        // Generate NRHS random solution vectors.

        var IX = 1;
        for (var J = 1; J <= NRHS; J++) {
          zlarnv(2, ISEED, N, XACT(IX));
          IX += LDA;
        }

        for (var ITRAN = 1; ITRAN <= 3; ITRAN++) {
          final TRANS = TRANSS[ITRAN - 1];
          final RCONDC = ITRAN == 1 ? RCONDO : RCONDI;

          // Set the right hand side.

          zlagtm(TRANS, N, NRHS, ONE, A, A(M + 1), A(N + M + 1),
              XACT.asMatrix(), LDA, ZERO, B.asMatrix(), LDA);

          // +    TEST 2
          // Solve op(A) * X = B and compute the residual.

          zlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);
          srnamc.SRNAMT = 'ZGTTRS';
          zgttrs(TRANS, N, NRHS, AF, AF(M + 1), AF(N + M + 1),
              AF(N + 2 * M + 1), IWORK, X.asMatrix(), LDA, INFO);

          // Check error code from ZGTTRS.

          if (INFO.value != 0) {
            alaerh(PATH, 'ZGTTRS', INFO.value, 0, TRANS, N, N, -1, -1, NRHS,
                IMAT, NFAIL, NERRS, NOUT);
          }

          zlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
          zgtt02(TRANS, N, NRHS, A, A(M + 1), A(N + M + 1), X.asMatrix(), LDA,
              WORK.asMatrix(), LDA, RESULT(2));

          // +    TEST 3
          // Check solution from generated exact solution.

          zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
              RESULT(3));

          // +    TESTS 4, 5, and 6
          // Use iterative refinement to improve the solution.

          srnamc.SRNAMT = 'ZGTRFS';
          zgtrfs(
              TRANS,
              N,
              NRHS,
              A,
              A(M + 1),
              A(N + M + 1),
              AF,
              AF(M + 1),
              AF(N + M + 1),
              AF(N + 2 * M + 1),
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

          // Check error code from ZGTRFS.

          if (INFO.value != 0) {
            alaerh(PATH, 'ZGTRFS', INFO.value, 0, TRANS, N, N, -1, -1, NRHS,
                IMAT, NFAIL, NERRS, NOUT);
          }

          zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
              RESULT(4));
          zgtt05(
              TRANS,
              N,
              NRHS,
              A,
              A(M + 1),
              A(N + M + 1),
              B.asMatrix(),
              LDA,
              X.asMatrix(),
              LDA,
              XACT.asMatrix(),
              LDA,
              RWORK,
              RWORK(NRHS + 1),
              RESULT(5));

          // Print information about the tests that did not pass the
          // threshold.

          for (var K = 2; K <= 6; K++) {
            if (RESULT[K] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.println(
                  ' TRANS=\'${TRANS.a1}\', N =${N.i5}, NRHS=${NRHS.i3}, type ${IMAT.i2}, test(${K.i2}) = ${RESULT[K].g12_5}');
              NFAIL++;
            }
          }
          NRUN += 5;
        }
      }
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
