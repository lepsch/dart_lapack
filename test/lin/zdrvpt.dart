import 'package:lapack/lapack.dart';

import '../matgen/zlatms.dart';
import 'aladhd.dart';
import 'alaerh.dart';
import 'alasvm.dart';
import 'common.dart';
import 'dget06.dart';
import 'zerrvx.dart';
import 'zget04.dart';
import 'zlaptm.dart';
import 'zlatb4.dart';
import 'zptt01.dart';
import 'zptt02.dart';
import 'zptt05.dart';

void zdrvpt(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NRHS,
  final double THRESH,
  final bool TSTERR,
  final Array<Complex> A_,
  final Array<double> D_,
  final Array<Complex> E_,
  final Array<Complex> B_,
  final Array<Complex> X_,
  final Array<Complex> XACT_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Nout NOUT,
) {
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

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0, ZERO = 0.0;
  const NTYPES = 12, NTESTS = 6;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS), Z = Array<double>(3);
  const ISEEDY = [0, 0, 0, 1];
  final INFO = Box(0);

  final PATH = '${'Zomplex precision'[0]}PT';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) zerrvx(PATH, NOUT);
  infoc.INFOT = 0;

  for (var IN = 1; IN <= NN; IN++) {
    // Do for each value of N in NVAL.

    final N = NVAL[IN];
    final LDA = max(1, N);
    final NIMAT = N <= 0 ? 1 : NTYPES;

    var IZERO = 0;

    for (var IMAT = 1; IMAT <= NIMAT; IMAT++) {
      // Do the tests only if DOTYPE( IMAT ) is true.

      if (N > 0 && !DOTYPE[IMAT]) continue;

      // Set up parameters with ZLATB4.

      final (:TYPE, :KL, :KU, :ANORM, :MODE, CNDNUM: COND, :DIST) =
          zlatb4(PATH, IMAT, N, N);

      final ZEROT = IMAT >= 8 && IMAT <= 10;
      if (IMAT <= 6) {
        // Type 1-6:  generate a symmetric tridiagonal matrix of
        // known condition number in lower triangular band storage.

        srnamc.SRNAMT = 'ZLATMS';
        zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, COND, ANORM, KL, KU, 'B',
            A.asMatrix(), 2, WORK, INFO);

        // Check the error code from ZLATMS.

        if (INFO.value != 0) {
          alaerh(PATH, 'ZLATMS', INFO.value, 0, ' ', N, N, KL, KU, -1, IMAT,
              NFAIL, NERRS, NOUT);
          continue;
        }
        IZERO = 0;

        // Copy the matrix to D and E.

        var IA = 1;
        for (var I = 1; I <= N - 1; I++) {
          D[I] = A[IA].real;
          E[I] = A[IA + 1];
          IA += 2;
        }
        if (N > 0) D[N] = A[IA].real;
      } else {
        // Type 7-12:  generate a diagonally dominant matrix with
        // unknown condition number in the vectors D and E.

        if (!ZEROT || !DOTYPE[7]) {
          // Let D and E have values from [-1,1].

          dlarnv(2, ISEED, N, D);
          zlarnv(2, ISEED, N - 1, E);

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
          if (N > 1) zdscal(N - 1, ANORM / DMAX, E, 1);
        } else if (IZERO > 0) {
          // Reuse the last matrix by copying back the zeroed out
          // elements.

          if (IZERO == 1) {
            D[1] = Z[2];
            if (N > 1) E[1] = Z[3].toComplex();
          } else if (IZERO == N) {
            E[N - 1] = Z[1].toComplex();
            D[N] = Z[2];
          } else {
            E[IZERO - 1] = Z[1].toComplex();
            D[IZERO] = Z[2];
            E[IZERO] = Z[3].toComplex();
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
            Z[3] = E[1].real;
            E[1] = Complex.zero;
          }
        } else if (IMAT == 9) {
          IZERO = N;
          if (N > 1) {
            Z[1] = E[N - 1].real;
            E[N - 1] = Complex.zero;
          }
          Z[2] = D[N];
          D[N] = ZERO;
        } else if (IMAT == 10) {
          IZERO = (N + 1) ~/ 2;
          if (IZERO > 1) {
            Z[1] = E[IZERO - 1].real;
            E[IZERO - 1] = Complex.zero;
            Z[3] = E[IZERO].real;
            E[IZERO] = Complex.zero;
          }
          Z[2] = D[IZERO];
          D[IZERO] = ZERO;
        }
      }

      // Generate NRHS random solution vectors.

      var IX = 1;
      for (var J = 1; J <= NRHS; J++) {
        zlarnv(2, ISEED, N, XACT(IX));
        IX += LDA;
      }

      // Set the right hand side.

      zlaptm('Lower', N, NRHS, ONE, D, E, XACT.asMatrix(), LDA, ZERO,
          B.asMatrix(), LDA);

      var RCONDC = ZERO;
      for (var IFACT = 1; IFACT <= 2; IFACT++) {
        final FACT = IFACT == 1 ? 'F' : 'N';

        // Compute the condition number for comparison with
        // the value returned by ZPTSVX.

        if (ZEROT) {
          if (IFACT == 1) continue;
          RCONDC = ZERO;
        } else if (IFACT == 1) {
          // Compute the 1-norm of A.

          final ANORM = zlanht('1', N, D, E);

          dcopy(N, D, 1, D(N + 1), 1);
          if (N > 1) zcopy(N - 1, E, 1, E(N + 1), 1);

          // Factor the matrix A.

          zpttrf(N, D(N + 1), E(N + 1), INFO);

          // Use ZPTTRS to solve for one column at a time of
          // inv(A), computing the maximum column sum as we go.

          var AINVNM = ZERO;
          for (var I = 1; I <= N; I++) {
            for (var J = 1; J <= N; J++) {
              X[J] = Complex.zero;
            }
            X[I] = Complex.one;
            zpttrs('Lower', N, 1, D(N + 1), E(N + 1), X.asMatrix(), LDA, INFO);
            AINVNM = max(AINVNM, dzasum(N, X, 1));
          }

          // Compute the 1-norm condition number of A.

          if (ANORM <= ZERO || AINVNM <= ZERO) {
            RCONDC = ONE;
          } else {
            RCONDC = (ONE / ANORM) / AINVNM;
          }
        }

        if (IFACT == 2) {
          // --- Test ZPTSV --

          dcopy(N, D, 1, D(N + 1), 1);
          if (N > 1) zcopy(N - 1, E, 1, E(N + 1), 1);
          zlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

          // Factor A as L*D*L' and solve the system A*X = B.

          srnamc.SRNAMT = 'ZPTSV';
          zptsv(N, NRHS, D(N + 1), E(N + 1), X.asMatrix(), LDA, INFO);

          // Check error code from ZPTSV .

          if (INFO.value != IZERO) {
            alaerh(PATH, 'ZPTSV ', INFO.value, IZERO, ' ', N, N, 1, 1, NRHS,
                IMAT, NFAIL, NERRS, NOUT);
          }
          final int NT;
          if (IZERO == 0) {
            // Check the factorization by computing the ratio
            //    norm(L*D*L' - A) / (n * norm(A) * EPS )

            zptt01(N, D, E, D(N + 1), E(N + 1), WORK, RESULT(1));

            // Compute the residual in the solution.

            zlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
            zptt02('Lower', N, NRHS, D, E, X.asMatrix(), LDA, WORK.asMatrix(),
                LDA, RESULT(2));

            // Check solution from generated exact solution.

            zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                RESULT(3));
            NT = 3;
          } else {
            NT = 0;
          }

          // Print information about the tests that did not pass
          // the threshold.

          for (var K = 1; K <= NT; K++) {
            if (RESULT[K] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
              NOUT.println(
                  ' ZPTSV , N =${N.i5}, type ${IMAT.i2}, test ${K.i2}, ratio = ${RESULT[K].g12_5}');
              NFAIL++;
            }
          }
          NRUN += NT;
        }

        // --- Test ZPTSVX ---

        if (IFACT > 1) {
          // Initialize D( N+1:2*N ) and E( N+1:2*N ) to zero.

          for (var I = 1; I <= N - 1; I++) {
            D[N + I] = ZERO;
            E[N + I] = Complex.zero;
          }
          if (N > 0) D[N + N] = ZERO;
        }

        zlaset('Full', N, NRHS, Complex.zero, Complex.zero, X.asMatrix(), LDA);

        // Solve the system and compute the condition number and
        // error bounds using ZPTSVX.

        final RCOND = Box(ZERO);
        srnamc.SRNAMT = 'ZPTSVX';
        zptsvx(
            FACT,
            N,
            NRHS,
            D,
            E,
            D(N + 1),
            E(N + 1),
            B.asMatrix(),
            LDA,
            X.asMatrix(),
            LDA,
            RCOND,
            RWORK,
            RWORK(NRHS + 1),
            WORK,
            RWORK(2 * NRHS + 1),
            INFO);

        // Check the error code from ZPTSVX.

        if (INFO.value != IZERO) {
          alaerh(PATH, 'ZPTSVX', INFO.value, IZERO, FACT, N, N, 1, 1, NRHS,
              IMAT, NFAIL, NERRS, NOUT);
        }
        final int K1;
        if (IZERO == 0) {
          if (IFACT == 2) {
            // Check the factorization by computing the ratio
            //    norm(L*D*L' - A) / (n * norm(A) * EPS )

            K1 = 1;
            zptt01(N, D, E, D(N + 1), E(N + 1), WORK, RESULT(1));
          } else {
            K1 = 2;
          }

          // Compute the residual in the solution.

          zlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
          zptt02('Lower', N, NRHS, D, E, X.asMatrix(), LDA, WORK.asMatrix(),
              LDA, RESULT(2));

          // Check solution from generated exact solution.

          zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
              RESULT(3));

          // Check error bounds from iterative refinement.

          zptt05(N, NRHS, D, E, B.asMatrix(), LDA, X.asMatrix(), LDA,
              XACT.asMatrix(), LDA, RWORK, RWORK(NRHS + 1), RESULT(4));
        } else {
          K1 = 6;
        }

        // Check the reciprocal of the condition number.

        RESULT[6] = dget06(RCOND.value, RCONDC);

        // Print information about the tests that did not pass
        // the threshold.

        for (var K = K1; K <= 6; K++) {
          if (RESULT[K] >= THRESH) {
            if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, PATH);
            NOUT.println(
                ' ZPTSVX, FACT=\'${FACT.a1}\', N =${N.i5}, type ${IMAT.i2}, test ${K.i2}, ratio = ${RESULT[K].g12_5}');
            NFAIL++;
          }
        }
        NRUN += 7 - K1;
      }
    }
  }

  // Print a summary of the results.

  alasvm(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
