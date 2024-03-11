import 'dart:math';

import 'package:lapack/src/blas/dasum.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlanst.dart';
import 'package:lapack/src/dlarnv.dart';
import 'package:lapack/src/dptcon.dart';
import 'package:lapack/src/dptrfs.dart';
import 'package:lapack/src/dpttrf.dart';
import 'package:lapack/src/dpttrs.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../matgen/dlatms.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'derrgt.dart';
import 'dget04.dart';
import 'dget06.dart';
import 'dlaptm.dart';
import 'dlatb4.dart';
import 'dptt01.dart';
import 'dptt02.dart';
import 'dptt05.dart';

void dchkpt(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NNS,
  final Array<int> NSVAL_,
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
) {
  final DOTYPE = DOTYPE_.having();
  final NVAL = NVAL_.having();
  final NSVAL = NSVAL_.having();
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
  const NTYPES = 12;
  const NTESTS = 7;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS), Z = Array<double>(3);
  const ISEEDY = [0, 0, 0, 1];
  final INFO = Box(0);

  final PATH = '${'Double precision'[0]}PT';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) derrgt(PATH, NOUT);
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

      // Set up parameters with DLATB4.

      final (:TYPE, :KL, :KU, ANORM: ANORMTMP, :MODE, :COND, :DIST) =
          dlatb4(PATH, IMAT, N, N);
      var ANORM = ANORMTMP;

      final ZEROT = IMAT >= 8 && IMAT <= 10;
      if (IMAT <= 6) {
        // Type 1-6:  generate a symmetric tridiagonal matrix of
        // known condition number in lower triangular band storage.

        srnamc.SRNAMT = 'DLATMS';
        dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, COND, ANORM, KL, KU, 'B',
            A.asMatrix(), 2, WORK, INFO);

        // Check the error code from DLATMS.

        if (INFO.value != 0) {
          alaerh(PATH, 'DLATMS', INFO.value, 0, ' ', N, N, KL, KU, -1, IMAT,
              NFAIL, NERRS, NOUT);
          continue;
        }
        IZERO = 0;

        // Copy the matrix to D and E.

        var IA = 1;
        for (var I = 1; I <= N - 1; I++) {
          D[I] = A[IA];
          E[I] = A[IA + 1];
          IA = IA + 2;
        }
        if (N > 0) D[N] = A[IA];
      } else {
        // Type 7-12:  generate a diagonally dominant matrix with
        // unknown condition number in the vectors D and E.

        if (!ZEROT || !DOTYPE[7]) {
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
          dscal(N - 1, ANORM / DMAX, E, 1);
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
            E[IZERO - 1] = ZERO;
            Z[3] = E[IZERO];
            E[IZERO] = ZERO;
          }
          Z[2] = D[IZERO];
          D[IZERO] = ZERO;
        }
      }

      dcopy(N, D, 1, D(N + 1), 1);
      if (N > 1) dcopy(N - 1, E, 1, E(N + 1), 1);

      // +    TEST 1
      // Factor A as L*D*L' and compute the ratio
      //    norm(L*D*L' - A) / (n * norm(A) * EPS )

      dpttrf(N, D(N + 1), E(N + 1), INFO);

      // Check error code from DPTTRF.

      if (INFO.value != IZERO) {
        alaerh(PATH, 'DPTTRF', INFO.value, IZERO, ' ', N, N, -1, -1, -1, IMAT,
            NFAIL, NERRS, NOUT);
        continue;
      }

      final double RCONDC;
      if (INFO.value > 0) {
        RCONDC = ZERO;
      } else {
        dptt01(N, D, E, D(N + 1), E(N + 1), WORK, RESULT(1));

        // Print the test ratio if greater than or equal to THRESH.

        if (RESULT[1] >= THRESH) {
          if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
          NOUT.println(
              ' N =${N.i5}, type ${IMAT.i2}, test ${1.i2}, ratio = ${RESULT[1].g12_5}');
          NFAIL++;
        }
        NRUN++;

        // Compute RCONDC = 1 / (norm(A) * norm(inv(A))

        // Compute norm(A).

        ANORM = dlanst('1', N, D, E);

        // Use DPTTRS to solve for one column at a time of inv(A),
        // computing the maximum column sum as we go.

        var AINVNM = ZERO;
        for (var I = 1; I <= N; I++) {
          for (var J = 1; J <= N; J++) {
            X[J] = ZERO;
          }
          X[I] = ONE;
          dpttrs(N, 1, D(N + 1), E(N + 1), X.asMatrix(), LDA, INFO);
          AINVNM = max(AINVNM, dasum(N, X, 1));
        }
        RCONDC = ONE / max(ONE, ANORM * AINVNM);

        for (var IRHS = 1; IRHS <= NNS; IRHS++) {
          final NRHS = NSVAL[IRHS];

          // Generate NRHS random solution vectors.

          var IX = 1;
          for (var J = 1; J <= NRHS; J++) {
            dlarnv(2, ISEED, N, XACT(IX));
            IX = IX + LDA;
          }

          // Set the right hand side.

          dlaptm(N, NRHS, ONE, D, E, XACT.asMatrix(), LDA, ZERO, B.asMatrix(),
              LDA);

          // +    TEST 2
          // Solve A*x = b and compute the residual.

          dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);
          dpttrs(N, NRHS, D(N + 1), E(N + 1), X.asMatrix(), LDA, INFO);

          // Check error code from DPTTRS.

          if (INFO.value != 0) {
            alaerh(PATH, 'DPTTRS', INFO.value, 0, ' ', N, N, -1, -1, NRHS, IMAT,
                NFAIL, NERRS, NOUT);
          }

          dlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
          dptt02(N, NRHS, D, E, X.asMatrix(), LDA, WORK.asMatrix(), LDA,
              RESULT(2));

          // +    TEST 3
          // Check solution from generated exact solution.

          dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
              RESULT(3));

          // +    TESTS 4, 5, and 6
          // Use iterative refinement to improve the solution.

          srnamc.SRNAMT = 'DPTRFS';
          dptrfs(N, NRHS, D, E, D(N + 1), E(N + 1), B.asMatrix(), LDA,
              X.asMatrix(), LDA, RWORK, RWORK(NRHS + 1), WORK, INFO);

          // Check error code from DPTRFS.

          if (INFO.value != 0) {
            alaerh(PATH, 'DPTRFS', INFO.value, 0, ' ', N, N, -1, -1, NRHS, IMAT,
                NFAIL, NERRS, NOUT);
          }

          dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
              RESULT(4));
          dptt05(N, NRHS, D, E, B.asMatrix(), LDA, X.asMatrix(), LDA,
              XACT.asMatrix(), LDA, RWORK, RWORK(NRHS + 1), RESULT(5));

          // Print information about the tests that did not pass the
          // threshold.

          for (var K = 2; K <= 6; K++) {
            if (RESULT[K] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.println(
                  ' N =${N.i5}, NRHS=${NRHS.i3}, type ${IMAT.i2}, test(${K.i2}) = ${(RESULT[K]).g12_5}');
              NFAIL++;
            }
          }
          NRUN += 5;
        }
      }

      // +    TEST 7
      // Estimate the reciprocal of the condition number of the
      // matrix.

      final RCOND = Box(0.0);
      srnamc.SRNAMT = 'DPTCON';
      dptcon(N, D(N + 1), E(N + 1), ANORM, RCOND, RWORK, INFO);

      // Check error code from DPTCON.

      if (INFO.value != 0) {
        alaerh(PATH, 'DPTCON', INFO.value, 0, ' ', N, N, -1, -1, -1, IMAT,
            NFAIL, NERRS, NOUT);
      }

      RESULT[7] = dget06(RCOND.value, RCONDC);

      // Print the test ratio if greater than or equal to THRESH.

      if (RESULT[7] >= THRESH) {
        if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
        NOUT.println(
            ' N =${N.i5}, type ${IMAT.i2}, test ${8.i2}, ratio = ${RESULT[7].g12_5}');
        NFAIL++;
      }
      NRUN++;
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
