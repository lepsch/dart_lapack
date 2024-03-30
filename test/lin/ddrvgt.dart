import 'dart:math';

import 'package:lapack/src/blas/dasum.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgtsv.dart';
import 'package:lapack/src/dgtsvx.dart';
import 'package:lapack/src/dgttrf.dart';
import 'package:lapack/src/dgttrs.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlagtm.dart';
import 'package:lapack/src/dlangt.dart';
import 'package:lapack/src/dlarnv.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../matgen/dlatms.dart';
import 'aladhd.dart';
import 'alaerh.dart';
import 'alasvm.dart';
import 'common.dart';
import 'derrvxx.dart';
import 'dget04.dart';
import 'dget06.dart';
import 'dgtt01.dart';
import 'dgtt02.dart';
import 'dgtt05.dart';
import 'dlatb4.dart';

void ddrvgt(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NRHS,
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
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DOTYPE = DOTYPE_.having();
  final NVAL = NVAL_.having();
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
  const NTESTS = 6;
  int NFAIL, NRUN;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS), Z = Array<double>(3);
  const ISEEDY = [0, 0, 0, 1], TRANSS = ['N', 'T', 'C'];
  final INFO = Box(0), NERRS = Box(0);

  final PATH = '${'Double precision'[0]}GT';
  NRUN = 0;
  NFAIL = 0;
  NERRS.value = 0;
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) derrvx(PATH, NOUT);
  infoc.INFOT = 0;

  for (var IN = 1; IN <= NN; IN++) {
    // Do for each value of N in NVAL.

    final N = NVAL[IN];
    final M = max(N - 1, 0);
    final LDA = max(1, N);
    final NIMAT = N <= 0 ? 1 : NTYPES;

    int IZERO = 0;
    for (var IMAT = 1; IMAT <= NIMAT; IMAT++) {
      // Do the tests only if DOTYPE( IMAT ) is true.

      if (!DOTYPE[IMAT]) continue;

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

        if (INFO.value != 0) {
          alaerh(PATH, 'DLATMS', INFO.value, 0, ' ', N, N, KL, KU, -1, IMAT,
              NFAIL, NERRS, NOUT);
          continue;
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

        if (!ZEROT || !DOTYPE[7]) {
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

      for (var IFACT = 1; IFACT <= 2; IFACT++) {
        final FACT = IFACT == 1 ? 'F' : 'N';

        // Compute the condition number for comparison with
        // the value returned by DGTSVX.
        double RCONDI = 0, RCONDO = 0;
        if (ZEROT) {
          if (IFACT == 1) continue;
          RCONDO = ZERO;
          RCONDI = ZERO;
        } else if (IFACT == 1) {
          dcopy(N + 2 * M, A, 1, AF, 1);

          // Compute the 1-norm and infinity-norm of A.

          final ANORMO = dlangt('1', N, A, A(M + 1), A(N + M + 1));
          final ANORMI = dlangt('I', N, A, A(M + 1), A(N + M + 1));

          // Factor the matrix A.

          dgttrf(
              N, AF, AF(M + 1), AF(N + M + 1), AF(N + 2 * M + 1), IWORK, INFO);

          // Use DGTTRS to solve for one column at a time of
          // inv(A), computing the maximum column sum as we go.

          var AINVNM = ZERO;
          for (var I = 1; I <= N; I++) {
            for (var J = 1; J <= N; J++) {
              X[J] = ZERO;
            }
            X[I] = ONE;
            dgttrs('No transpose', N, 1, AF, AF(M + 1), AF(N + M + 1),
                AF(N + 2 * M + 1), IWORK, X.asMatrix(), LDA, INFO);
            AINVNM = max(AINVNM, dasum(N, X, 1));
          }

          // Compute the 1-norm condition number of A.

          if (ANORMO <= ZERO || AINVNM <= ZERO) {
            RCONDO = ONE;
          } else {
            RCONDO = (ONE / ANORMO) / AINVNM;
          }

          // Use DGTTRS to solve for one column at a time of
          // inv(A'), computing the maximum column sum as we go.

          AINVNM = ZERO;
          for (var I = 1; I <= N; I++) {
            for (var J = 1; J <= N; J++) {
              X[J] = ZERO;
            }
            X[I] = ONE;
            dgttrs('Transpose', N, 1, AF, AF(M + 1), AF(N + M + 1),
                AF(N + 2 * M + 1), IWORK, X.asMatrix(), LDA, INFO);
            AINVNM = max(AINVNM, dasum(N, X, 1));
          }

          // Compute the infinity-norm condition number of A.

          if (ANORMI <= ZERO || AINVNM <= ZERO) {
            RCONDI = ONE;
          } else {
            RCONDI = (ONE / ANORMI) / AINVNM;
          }
        }

        int NT = 0;
        for (var ITRAN = 1; ITRAN <= 3; ITRAN++) {
          final TRANS = TRANSS[ITRAN - 1];
          final RCONDC = ITRAN == 1 ? RCONDO : RCONDI;

          // Generate NRHS random solution vectors.

          var IX = 1;
          for (var J = 1; J <= NRHS; J++) {
            dlarnv(2, ISEED, N, XACT(IX));
            IX += LDA;
          }

          // Set the right hand side.

          dlagtm(TRANS, N, NRHS, ONE, A, A(M + 1), A(N + M + 1),
              XACT.asMatrix(), LDA, ZERO, B.asMatrix(), LDA);

          if (IFACT == 2 && ITRAN == 1) {
            // --- Test DGTSV  ---

            // Solve the system using Gaussian elimination with
            // partial pivoting.

            dcopy(N + 2 * M, A, 1, AF, 1);
            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

            srnamc.SRNAMT = 'DGTSV ';
            dgtsv(
                N, NRHS, AF, AF(M + 1), AF(N + M + 1), X.asMatrix(), LDA, INFO);

            // Check error code from DGTSV .

            if (INFO.value != IZERO) {
              alaerh(PATH, 'DGTSV ', INFO.value, IZERO, ' ', N, N, 1, 1, NRHS,
                  IMAT, NFAIL, NERRS, NOUT);
            }
            NT = 1;
            if (IZERO == 0) {
              // Check residual of computed solution.

              dlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
              dgtt02(TRANS, N, NRHS, A, A(M + 1), A(N + M + 1), X.asMatrix(),
                  LDA, WORK.asMatrix(), LDA, RESULT(2));

              // Check solution from generated exact solution.

              dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                  RESULT(3));
              NT = 3;
            }

            // Print information about the tests that did not pass
            // the threshold.

            for (var K = 2; K <= NT; K++) {
              if (RESULT[K] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
                NOUT.println(
                    ' DGTSV , N =${N.i5}, type ${IMAT.i2}, test ${K.i2}, ratio = ${RESULT[K].g12_5}');
                NFAIL++;
              }
            }
            NRUN += NT - 1;
          }

          // --- Test DGTSVX ---

          if (IFACT > 1) {
            // Initialize AF to zero.

            for (var I = 1; I <= 3 * N - 2; I++) {
              AF[I] = ZERO;
            }
          }
          dlaset('Full', N, NRHS, ZERO, ZERO, X.asMatrix(), LDA);

          // Solve the system and compute the condition number and
          // error bounds using DGTSVX.

          srnamc.SRNAMT = 'DGTSVX';
          final RCOND = Box(0.0);
          dgtsvx(
              FACT,
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
              RCOND,
              RWORK,
              RWORK(NRHS + 1),
              WORK,
              IWORK(N + 1),
              INFO);

          // Check the error code from DGTSVX.

          if (INFO.value != IZERO) {
            alaerh(PATH, 'DGTSVX', INFO.value, IZERO, FACT + TRANS, N, N, 1, 1,
                NRHS, IMAT, NFAIL, NERRS, NOUT);
          }

          final int K1;
          if (IFACT >= 2) {
            // Reconstruct matrix from factors and compute
            // residual.

            dgtt01(
                N,
                A,
                A(M + 1),
                A(N + M + 1),
                AF,
                AF(M + 1),
                AF(N + M + 1),
                AF(N + 2 * M + 1),
                IWORK,
                WORK.asMatrix(),
                LDA,
                RWORK,
                RESULT(1));
            K1 = 1;
          } else {
            K1 = 2;
          }

          if (INFO.value == 0) {
            // Check residual of computed solution.

            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
            dgtt02(TRANS, N, NRHS, A, A(M + 1), A(N + M + 1), X.asMatrix(), LDA,
                WORK.asMatrix(), LDA, RESULT(2));

            // Check solution from generated exact solution.

            dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                RESULT(3));

            // Check the error bounds from iterative refinement.

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
                RESULT(4));
            NT = 5;
          }

          // Print information about the tests that did not pass
          // the threshold.

          int K;
          for (K = K1; K <= NT; K++) {
            if (RESULT[K] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
              NOUT.print9998('DGTSVX', FACT, TRANS, N, IMAT, K, RESULT[K]);
              NFAIL++;
            }
          }

          // Check the reciprocal of the condition number.

          RESULT[6] = dget06(RCOND.value, RCONDC);
          if (RESULT[6] >= THRESH) {
            if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
            NOUT.print9998('DGTSVX', FACT, TRANS, N, IMAT, K, RESULT[K]);
            NFAIL++;
          }
          NRUN += NT - K1 + 2;
        }
      }
    }
  }

  // Print a summary of the results.

  alasvm(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}

extension on Nout {
  void print9998(String s, String fact, String trans, int n, int type, int test,
      double ratio) {
    println(
        ' $s, FACT=\'${fact.a1}\', TRANS=\'${trans.a1}\', N =${n.i5}, type ${type.i2}, test ${test.i2}, ratio = ${ratio.g12_5}');
  }
}
