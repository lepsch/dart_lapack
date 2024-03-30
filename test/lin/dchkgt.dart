import 'dart:math';

import 'package:lapack/lapack.dart';
import 'package:test/test.dart';

import '../matgen/dlatms.dart';
import '../test_driver.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'derrge.dart';
import 'dget04.dart';
import 'dget06.dart';
import 'dgtt01.dart';
import 'dgtt02.dart';
import 'dgtt05.dart';
import 'dlatb4.dart';

void dchkgt(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NNS,
  final Array<int> NSVAL_,
  final double THRESH,
  final bool TSTERR,
  final Array<double> A_,
  final Array<double> AF_,
  final Array<double> B_,
  final Array<double> X_,
  final Array<double> XACT_,
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
  const NTYPES = 12, NTESTS = 7;
  final RESULT = Array<double>(NTESTS), Z = Array<double>(3);
  const ISEEDY = [0, 0, 0, 1], TRANSS = ['N', 'T', 'C'];

  final PATH = '${'Double precision'[0]}GT';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  final ISEED = Array.fromList(ISEEDY);

  test.group('error exits', () {
    // Test the error exits
    if (TSTERR) derrge(PATH, NOUT, test);
    test.tearDown(() {
      infoc.INFOT = 0;
    });
  });

  for (final IN in 1.through(NN)) {
    // Do for each value of N in NVAL.

    final N = NVAL[IN];
    final M = max(N - 1, 0);
    final LDA = max(1, N);
    final NIMAT = N <= 0 ? 1 : NTYPES;

    int IZERO = 0;
    for (final IMAT in 1.through(NIMAT)) {
      // Do the tests only if DOTYPE( IMAT ) is true.
      final skip = !DOTYPE[IMAT];

      test('DCHKGT (IN=$IN IMAT=$IMAT)', () {
        final INFO = Box(0);

        // Set up parameters with DLATB4.

        final (:TYPE, :KL, :KU, :ANORM, :MODE, :COND, :DIST) =
            dlatb4(PATH, IMAT, N, N);

        final ZEROT = IMAT >= 8 && IMAT <= 10;
        if (IMAT <= 6) {
          // Types 1-6:  generate matrices of known condition number.

          final KOFF = max(2 - KU, 3 - max(1, N)).toInt();
          srnamc.SRNAMT = 'DLATMS';
          dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, COND, ANORM, KL, KU, 'Z',
              AF(KOFF).asMatrix(), 3, WORK, INFO);

          // Check the error code from DLATMS.
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            alaerh(PATH, 'DLATMS', INFO.value, 0, ' ', N, N, KL, KU, -1, IMAT,
                NFAIL, NERRS, NOUT);
            return;
          }
          IZERO = 0;

          if (N > 1) {
            dcopy(N - 1, AF(4), 3, A, 1);
            dcopy(N - 1, AF(3), 3, A(N + M + 1), 1);
          }
          dcopy(N, AF(2), 3, A(M + 1), 1);
        } else {
          // Types 7-12:  generate tridiagonal matrices with
          // unknown condition numbers.

          if (!ZEROT || !DOTYPE[7] || TestDriver.isAsync) {
            // Generate a matrix with elements from [-1,1].

            dlarnv(2, ISEED, N + 2 * M, A);
            if (ANORM != ONE) dscal(N + 2 * M, ANORM, A, 1);
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
            A[N] = ZERO;
            if (N > 1) {
              Z[3] = A[1];
              A[1] = ZERO;
            }
          } else if (IMAT == 9) {
            IZERO = N;
            Z[1] = A[3 * N - 2];
            Z[2] = A[2 * N - 1];
            A[3 * N - 2] = ZERO;
            A[2 * N - 1] = ZERO;
          } else {
            IZERO = (N + 1) ~/ 2;
            for (var I = IZERO; I <= N - 1; I++) {
              A[2 * N - 2 + I] = ZERO;
              A[N - 1 + I] = ZERO;
              A[I] = ZERO;
            }
            A[3 * N - 2] = ZERO;
            A[2 * N - 1] = ZERO;
          }
        }

        // +    TEST 1
        // Factor A as L*U and compute the ratio
        //    norm(L*U - A) / (n * norm(A) * EPS )

        dcopy(N + 2 * M, A, 1, AF, 1);
        srnamc.SRNAMT = 'DGTTRF';
        dgttrf(N, AF, AF(M + 1), AF(N + M + 1), AF(N + 2 * M + 1), IWORK, INFO);

        // Check error code from DGTTRF.
        test.expect(INFO.value, IZERO);
        if (INFO.value != IZERO) {
          alaerh(PATH, 'DGTTRF', INFO.value, IZERO, ' ', N, N, 1, 1, -1, IMAT,
              NFAIL, NERRS, NOUT);
        }
        final TRFCON = INFO.value != 0;

        dgtt01(N, A, A(M + 1), A(N + M + 1), AF, AF(M + 1), AF(N + M + 1),
            AF(N + 2 * M + 1), IWORK, WORK.asMatrix(), LDA, RWORK, RESULT(1));

        // Print the test ratio if it is >= THRESH.
        final reason =
            '${' ' * 12}N =${N.i5},${' ' * 10} type ${IMAT.i2}, test(${1.i2}) = ${RESULT[1].g12_5}';
        test.expect(RESULT[1], lessThan(THRESH), reason: reason);
        if (RESULT[1] >= THRESH) {
          if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
          NOUT.println(reason);
          NFAIL++;
        }
        NRUN++;

        double RCONDI = ZERO, RCONDO = ZERO;
        for (var ITRAN = 1; ITRAN <= 2; ITRAN++) {
          final TRANS = TRANSS[ITRAN - 1];
          final NORM = ITRAN == 1 ? 'O' : 'I';

          final ANORM = dlangt(NORM, N, A, A(M + 1), A(N + M + 1));

          final double RCONDC;
          if (!TRFCON) {
            // Use DGTTRS to solve for one column at a time of inv(A)
            // or inv(A^T), computing the maximum column sum as we
            // go.

            var AINVNM = ZERO;
            for (var I = 1; I <= N; I++) {
              for (var J = 1; J <= N; J++) {
                X[J] = ZERO;
              }
              X[I] = ONE;
              dgttrs(TRANS, N, 1, AF, AF(M + 1), AF(N + M + 1),
                  AF(N + 2 * M + 1), IWORK, X.asMatrix(), LDA, INFO);
              AINVNM = max(AINVNM, dasum(N, X, 1));
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

          srnamc.SRNAMT = 'DGTCON';
          final RCOND = Box(0.0);
          dgtcon(NORM, N, AF, AF(M + 1), AF(N + M + 1), AF(N + 2 * M + 1),
              IWORK, ANORM, RCOND, WORK, IWORK(N + 1), INFO);

          // Check error code from DGTCON.
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            alaerh(PATH, 'DGTCON', INFO.value, 0, NORM, N, N, -1, -1, -1, IMAT,
                NFAIL, NERRS, NOUT);
          }

          RESULT[7] = dget06(RCOND.value, RCONDC);

          // Print the test ratio if it is >= THRESH.
          final reason =
              ' NORM =\'${NORM.a1}\', N =${N.i5},${' ' * 10} type ${IMAT.i2}, test(${7.i2}) = ${RESULT[7].g12_5}';
          test.expect(RESULT[7], lessThan(THRESH), reason: reason);
          if (RESULT[7] >= THRESH) {
            if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
            NOUT.println(reason);
            NFAIL++;
          }
          NRUN++;
        }

        // Skip the remaining tests if the matrix is singular.
        if (TRFCON) return;

        for (var IRHS = 1; IRHS <= NNS; IRHS++) {
          final NRHS = NSVAL[IRHS];

          // Generate NRHS random solution vectors.

          var IX = 1;
          for (var J = 1; J <= NRHS; J++) {
            dlarnv(2, ISEED, N, XACT(IX));
            IX += LDA;
          }

          for (var ITRAN = 1; ITRAN <= 3; ITRAN++) {
            final TRANS = TRANSS[ITRAN - 1];
            final RCONDC = ITRAN == 1 ? RCONDO : RCONDI;

            // Set the right hand side.

            dlagtm(TRANS, N, NRHS, ONE, A, A(M + 1), A(N + M + 1),
                XACT.asMatrix(), LDA, ZERO, B.asMatrix(), LDA);

            // +    TEST 2
            // Solve op(A) * X = B and compute the residual.

            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);
            srnamc.SRNAMT = 'DGTTRS';
            dgttrs(TRANS, N, NRHS, AF, AF(M + 1), AF(N + M + 1),
                AF(N + 2 * M + 1), IWORK, X.asMatrix(), LDA, INFO);

            // Check error code from DGTTRS.
            test.expect(INFO.value, 0);
            if (INFO.value != 0) {
              alaerh(PATH, 'DGTTRS', INFO.value, 0, TRANS, N, N, -1, -1, NRHS,
                  IMAT, NFAIL, NERRS, NOUT);
            }

            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
            dgtt02(TRANS, N, NRHS, A, A(M + 1), A(N + M + 1), X.asMatrix(), LDA,
                WORK.asMatrix(), LDA, RESULT(2));

            // +    TEST 3
            // Check solution from generated exact solution.

            dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                RESULT(3));

            // +    TESTS 4, 5, and 6
            // Use iterative refinement to improve the solution.

            srnamc.SRNAMT = 'DGTRFS';
            dgtrfs(
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
                IWORK(N + 1),
                INFO);

            // Check error code from DGTRFS.
            test.expect(INFO.value, 0);
            if (INFO.value != 0) {
              alaerh(PATH, 'DGTRFS', INFO.value, 0, TRANS, N, N, -1, -1, NRHS,
                  IMAT, NFAIL, NERRS, NOUT);
            }

            dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                RESULT(4));
            dgtt05(
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

            // Print information about the tests that did not pass
            // the threshold.

            for (var K = 2; K <= 6; K++) {
              final reason =
                  ' TRANS=\'${TRANS.a1}\', N =${N.i5}, NRHS=${NRHS.i3}, type ${IMAT.i2}, test(${K.i2}) = ${RESULT[K].g12_5}';
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
      }, skip: skip);
    }
  }

  // Print a summary of the results.
  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
