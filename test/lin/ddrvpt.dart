import 'dart:math';

import 'package:lapack/lapack.dart';
import 'package:test/test.dart';

import '../matgen/dlatms.dart';
import '../test_driver.dart';
import 'aladhd.dart';
import 'alaerh.dart';
import 'alasvm.dart';
import 'common.dart';
import 'derrvx.dart';
import 'dget04.dart';
import 'dget06.dart';
import 'dlaptm.dart';
import 'dlatb4.dart';
import 'dptt01.dart';
import 'dptt02.dart';
import 'dptt05.dart';

void ddrvpt(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NRHS,
  final double THRESH,
  final bool TSTERR,
  final Array<double> A_,
  final Array<double> D_,
  final Array<double> E_,
  final Array<double> B_,
  final Array<double> X_,
  final Array<double> XACT_,
  final Array<double> WORK_,
  final Array<double> RWORK_,
  final Nout NOUT,
  final TestDriver test,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DOTYPE = DOTYPE_.having();
  final NVAL = NVAL_.having();
  final A = A_.having();
  final D = D_.having();
  final E = E_.having();
  final B = B_.having();
  final X = X_.having();
  final XACT = XACT_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  const NTYPES = 12, NTESTS = 6;
  final RESULT = Array<double>(NTESTS), Z = Array<double>(3);
  const ISEEDY = [0, 0, 0, 1];

  final PATH = '${'Double precision'[0]}PT';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  final ISEED = Array.fromList(ISEEDY);

  test.group('error exits', () {
    // Test the error exits
    if (TSTERR) derrvx(PATH, NOUT, test);
    test.tearDown(() {
      infoc.INFOT = 0;
    });
  });

  for (final IN in 1.through(NN)) {
    // Do for each value of N in NVAL.

    final N = NVAL[IN];
    final LDA = max(1, N);
    final NIMAT = N <= 0 ? 1 : NTYPES;
    var IZERO = 0;

    for (final IMAT in 1.through(NIMAT)) {
      // Do the tests only if DOTYPE( IMAT ) is true.
      final skip = N > 0 && !DOTYPE[IMAT];

      test('DDRVPT (IN=$IN IMAT=$IMAT)', () {
        final INFO = Box(0);

        // Set up parameters with DLATB4.
        final (:TYPE, :KL, :KU, :ANORM, :MODE, :COND, :DIST) =
            dlatb4(PATH, IMAT, N, N);

        final ZEROT = IMAT >= 8 && IMAT <= 10;
        if (IMAT <= 6) {
          // Type 1-6:  generate a symmetric tridiagonal matrix of
          // known condition number in lower triangular band storage.

          srnamc.SRNAMT = 'DLATMS';
          dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, COND, ANORM, KL, KU, 'B',
              A.asMatrix(), 2, WORK, INFO);

          // Check the error code from DLATMS.
          test.expect(INFO.value, 0);
          if (INFO.value != 0) {
            alaerh(PATH, 'DLATMS', INFO.value, 0, ' ', N, N, KL, KU, -1, IMAT,
                NFAIL, NERRS, NOUT);
            return;
          }
          IZERO = 0;

          // Copy the matrix to D and E.

          var IA = 1;
          for (var I = 1; I <= N - 1; I++) {
            D[I] = A[IA];
            E[I] = A[IA + 1];
            IA += 2;
          }
          if (N > 0) D[N] = A[IA];
        } else {
          // Type 7-12:  generate a diagonally dominant matrix with
          // unknown condition number in the vectors D and E.

          if (!ZEROT || !DOTYPE[7] || TestDriver.isAsync) {
            // Let D and E have values from [-1,1].

            dlarnv(2, ISEED, N, D);
            dlarnv(2, ISEED, N - 1, E);

            // Make the tridiagonal matrix diagonally dominant.

            if (N == 1) {
              D[1] = D[1].abs();
            } else {
              D[1] = D[1].abs() + E[1].abs();
              D[N] = D[N].abs() + E[N - 1].abs();
              for (var I = 2; I <= N - 1; I++) {
                D[I] = D[I].abs() + E[I].abs() + E[I - 1].abs();
              }
            }

            // Scale D and E so the maximum element is ANORM.

            final IX = idamax(N, D, 1);
            final DMAX = D[IX];
            dscal(N, ANORM / DMAX, D, 1);
            if (N > 1) dscal(N - 1, ANORM / DMAX, E, 1);
          } else if (IZERO > 0) {
            // Reuse the last matrix by copying back the zeroed out
            // elements.

            if (IZERO == 1) {
              D[1] = Z[2];
              if (N > 1) E[1] = Z[3];
            } else if (IZERO == N) {
              E[N - 1] = Z[1];
              D[N] = Z[2];
            } else {
              E[IZERO - 1] = Z[1];
              D[IZERO] = Z[2];
              E[IZERO] = Z[3];
            }
          }

          // For types 8-10, set one row and column of the matrix to
          // zero.

          IZERO = 0;
          if (IMAT == 8) {
            IZERO = 1;
            Z[2] = D[1];
            D[1] = ZERO;
            if (N > 1) {
              Z[3] = E[1];
              E[1] = ZERO;
            }
          } else if (IMAT == 9) {
            IZERO = N;
            if (N > 1) {
              Z[1] = E[N - 1];
              E[N - 1] = ZERO;
            }
            Z[2] = D[N];
            D[N] = ZERO;
          } else if (IMAT == 10) {
            IZERO = (N + 1) ~/ 2;
            if (IZERO > 1) {
              Z[1] = E[IZERO - 1];
              Z[3] = E[IZERO];
              E[IZERO - 1] = ZERO;
              E[IZERO] = ZERO;
            }
            Z[2] = D[IZERO];
            D[IZERO] = ZERO;
          }
        }

        // Generate NRHS random solution vectors.

        var IX = 1;
        for (var J = 1; J <= NRHS; J++) {
          dlarnv(2, ISEED, N, XACT(IX));
          IX += LDA;
        }

        // Set the right hand side.

        dlaptm(
            N, NRHS, ONE, D, E, XACT.asMatrix(), LDA, ZERO, B.asMatrix(), LDA);
        var RCONDC = ZERO;
        for (var IFACT = 1; IFACT <= 2; IFACT++) {
          final FACT = IFACT == 1 ? 'F' : 'N';

          // Compute the condition number for comparison with
          // the value returned by DPTSVX.

          if (ZEROT) {
            if (IFACT == 1) continue;
            RCONDC = ZERO;
          } else if (IFACT == 1) {
            // Compute the 1-norm of A.

            final ANORM = dlanst('1', N, D, E);

            dcopy(N, D, 1, D(N + 1), 1);
            if (N > 1) dcopy(N - 1, E, 1, E(N + 1), 1);

            // Factor the matrix A.

            dpttrf(N, D(N + 1), E(N + 1), INFO);

            // Use DPTTRS to solve for one column at a time of
            // inv(A), computing the maximum column sum as we go.

            var AINVNM = ZERO;
            for (var I = 1; I <= N; I++) {
              for (var J = 1; J <= N; J++) {
                X[J] = ZERO;
              }
              X[I] = ONE;
              dpttrs(N, 1, D(N + 1), E(N + 1), X.asMatrix(), LDA, INFO);
              AINVNM = max(AINVNM, dasum(N, X, 1));
            }

            // Compute the 1-norm condition number of A.

            if (ANORM <= ZERO || AINVNM <= ZERO) {
              RCONDC = ONE;
            } else {
              RCONDC = (ONE / ANORM) / AINVNM;
            }
          }

          if (IFACT == 2) {
            // --- Test DPTSV --

            dcopy(N, D, 1, D(N + 1), 1);
            if (N > 1) dcopy(N - 1, E, 1, E(N + 1), 1);
            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

            // Factor A as L*D*L' and solve the system A*X = B.

            srnamc.SRNAMT = 'DPTSV ';
            dptsv(N, NRHS, D(N + 1), E(N + 1), X.asMatrix(), LDA, INFO);

            // Check error code from DPTSV .
            test.expect(INFO.value, IZERO);
            if (INFO.value != IZERO) {
              alaerh(PATH, 'DPTSV ', INFO.value, IZERO, ' ', N, N, 1, 1, NRHS,
                  IMAT, NFAIL, NERRS, NOUT);
            }
            final int NT;
            if (IZERO == 0) {
              // Check the factorization by computing the ratio
              //    norm(L*D*L' - A) / (n * norm(A) * EPS )

              dptt01(N, D, E, D(N + 1), E(N + 1), WORK, RESULT(1));

              // Compute the residual in the solution.

              dlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
              dptt02(N, NRHS, D, E, X.asMatrix(), LDA, WORK.asMatrix(), LDA,
                  RESULT(2));

              // Check solution from generated exact solution.

              dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                  RESULT(3));
              NT = 3;
            } else {
              NT = 0;
            }

            // Print information about the tests that did not pass
            // the threshold.

            for (var K = 1; K <= NT; K++) {
              final reason =
                  ' DPTSV, N =${N.i5}, type ${IMAT.i2}, test ${K.i2}, ratio = ${RESULT[K].g12_5}';
              test.expect(RESULT[K], lessThan(THRESH), reason: reason);
              if (RESULT[K] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                NOUT.println(reason);
                NFAIL++;
              }
            }
            NRUN += NT;
          }

          // --- Test DPTSVX ---

          if (IFACT > 1) {
            // Initialize D( N+1:2*N ) and E( N+1:2*N ) to zero.

            for (var I = 1; I <= N - 1; I++) {
              D[N + I] = ZERO;
              E[N + I] = ZERO;
            }
            if (N > 0) D[N + N] = ZERO;
          }

          dlaset('Full', N, NRHS, ZERO, ZERO, X.asMatrix(), LDA);

          // Solve the system and compute the condition number and
          // error bounds using DPTSVX.

          final RCOND = Box(0.0);
          srnamc.SRNAMT = 'DPTSVX';
          dptsvx(FACT, N, NRHS, D, E, D(N + 1), E(N + 1), B.asMatrix(), LDA,
              X.asMatrix(), LDA, RCOND, RWORK, RWORK(NRHS + 1), WORK, INFO);

          // Check the error code from DPTSVX.
          test.expect(INFO.value, IZERO);
          if (INFO.value != IZERO) {
            alaerh(PATH, 'DPTSVX', INFO.value, IZERO, FACT, N, N, 1, 1, NRHS,
                IMAT, NFAIL, NERRS, NOUT);
          }
          final int K1;
          if (IZERO == 0) {
            if (IFACT == 2) {
              // Check the factorization by computing the ratio
              //    norm(L*D*L' - A) / (n * norm(A) * EPS )

              K1 = 1;
              dptt01(N, D, E, D(N + 1), E(N + 1), WORK, RESULT(1));
            } else {
              K1 = 2;
            }

            // Compute the residual in the solution.

            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
            dptt02(N, NRHS, D, E, X.asMatrix(), LDA, WORK.asMatrix(), LDA,
                RESULT(2));

            // Check solution from generated exact solution.

            dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                RESULT(3));

            // Check error bounds from iterative refinement.

            dptt05(N, NRHS, D, E, B.asMatrix(), LDA, X.asMatrix(), LDA,
                XACT.asMatrix(), LDA, RWORK, RWORK(NRHS + 1), RESULT(4));
          } else {
            K1 = 6;
          }

          // Check the reciprocal of the condition number.

          RESULT[6] = dget06(RCOND.value, RCONDC);

          // Print information about the tests that did not pass
          // the threshold.

          for (var K = K1; K <= 6; K++) {
            final reason =
                ' DPTSVX, FACT=\'${FACT.a1}\', N =${N.i5}, type ${IMAT.i2}, test ${K.i2}, ratio = ${RESULT[K].g12_5}';
            test.expect(RESULT[K], lessThan(THRESH), reason: reason);
            if (RESULT[K] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
              NOUT.println(reason);
              NFAIL++;
            }
          }
          NRUN += 7 - K1;
        }
      }, skip: skip);
    }
  }

  // Print a summary of the results.
  alasvm(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
