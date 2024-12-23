// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/blas/dzasum.dart';
import 'package:dart_lapack/src/blas/idamax.dart';
import 'package:dart_lapack/src/blas/zcopy.dart';
import 'package:dart_lapack/src/blas/zdscal.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/dlarnv.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlanht.dart';
import 'package:dart_lapack/src/zlarnv.dart';
import 'package:dart_lapack/src/zptcon.dart';
import 'package:dart_lapack/src/zptrfs.dart';
import 'package:dart_lapack/src/zpttrf.dart';
import 'package:dart_lapack/src/zpttrs.dart';

import '../matgen/zlatms.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'dget06.dart';
import 'zerrgt.dart';
import 'zget04.dart';
import 'zlaptm.dart';
import 'zlatb4.dart';
import 'zptt01.dart';
import 'zptt02.dart';
import 'zptt05.dart';

void zchkpt(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NNS,
  final Array<int> NSVAL_,
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
  final NSVAL = NSVAL_.having();
  final A = A_.having();
  final D = D_.having();
  final E = E_.having();
  final B = B_.having();
  final X = X_.having();
  final XACT = XACT_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  const NTYPES = 12, NTESTS = 7;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  final Z = Array<Complex>(3);
  const ISEEDY = [0, 0, 0, 1], UPLOS = ['U', 'L'];
  final INFO = Box(0);

  final PATH = '${'Zomplex precision'[0]}PT';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) zerrgt(PATH, NOUT);
  infoc.INFOT = 0;

  for (var IN = 1; IN <= NN; IN++) {
    // Do for each value of N in NVAL.

    final N = NVAL[IN];
    final LDA = max(1, N);
    final NIMAT = N <= 0 ? 1 : NTYPES;

    var IZERO = 0, RCONDC = ZERO;
    for (var IMAT = 1; IMAT <= NIMAT; IMAT++) {
      // Do the tests only if DOTYPE( IMAT ) is true.

      if (N > 0 && !DOTYPE[IMAT]) continue;

      // Set up parameters with ZLATB4.

      final (:TYPE, :KL, :KU, :ANORM, :MODE, CNDNUM: COND, :DIST) =
          zlatb4(PATH, IMAT, N, N);

      final ZEROT = IMAT >= 8 && IMAT <= 10;
      if (IMAT <= 6) {
        // Type 1-6:  generate a Hermitian tridiagonal matrix of
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
          // Let E be complex, D real, with values from [-1,1].

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
          zdscal(N - 1, ANORM / DMAX, E, 1);
        } else if (IZERO > 0) {
          // Reuse the last matrix by copying back the zeroed out
          // elements.

          if (IZERO == 1) {
            D[1] = Z[2].real;
            if (N > 1) E[1] = Z[3];
          } else if (IZERO == N) {
            E[N - 1] = Z[1];
            D[N] = Z[2].real;
          } else {
            E[IZERO - 1] = Z[1];
            D[IZERO] = Z[2].real;
            E[IZERO] = Z[3];
          }
        }

        // For types 8-10, set one row and column of the matrix to
        // zero.

        IZERO = 0;
        if (IMAT == 8) {
          IZERO = 1;
          Z[2] = D[1].toComplex();
          D[1] = ZERO;
          if (N > 1) {
            Z[3] = E[1];
            E[1] = Complex.zero;
          }
        } else if (IMAT == 9) {
          IZERO = N;
          if (N > 1) {
            Z[1] = E[N - 1];
            E[N - 1] = Complex.zero;
          }
          Z[2] = D[N].toComplex();
          D[N] = ZERO;
        } else if (IMAT == 10) {
          IZERO = (N + 1) ~/ 2;
          if (IZERO > 1) {
            Z[1] = E[IZERO - 1];
            Z[3] = E[IZERO];
            E[IZERO - 1] = Complex.zero;
            E[IZERO] = Complex.zero;
          }
          Z[2] = D[IZERO].toComplex();
          D[IZERO] = ZERO;
        }
      }

      dcopy(N, D, 1, D(N + 1), 1);
      if (N > 1) zcopy(N - 1, E, 1, E(N + 1), 1);

      // +    TEST 1
      // Factor A as L*D*L' and compute the ratio
      //    norm(L*D*L' - A) / (n * norm(A) * EPS )

      zpttrf(N, D(N + 1), E(N + 1), INFO);

      // Check error code from ZPTTRF.

      if (INFO.value != IZERO) {
        alaerh(PATH, 'ZPTTRF', INFO.value, IZERO, ' ', N, N, -1, -1, -1, IMAT,
            NFAIL, NERRS, NOUT);
        continue;
      }

      if (INFO.value > 0) {
        RCONDC = ZERO;
      } else {
        zptt01(N, D, E, D(N + 1), E(N + 1), WORK, RESULT(1));

        // Print the test ratio if greater than or equal to THRESH.

        if (RESULT[1] >= THRESH) {
          if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
          NOUT.print9999(N, IMAT, 1, RESULT[1]);
          NFAIL++;
        }
        NRUN++;

        // Compute RCONDC = 1 / (norm(A) * norm(inv(A))

        // Compute norm(A).

        final ANORM = zlanht('1', N, D, E);

        // Use ZPTTRS to solve for one column at a time of inv(A),
        // computing the maximum column sum as we go.

        var AINVNM = ZERO;
        for (var I = 1; I <= N; I++) {
          for (var J = 1; J <= N; J++) {
            X[J] = Complex.zero;
          }
          X[I] = Complex.one;
          zpttrs('Lower', N, 1, D(N + 1), E(N + 1), X.asMatrix(), LDA, INFO);
          AINVNM = max(AINVNM, dzasum(N, X, 1));
        }
        RCONDC = ONE / max(ONE, ANORM * AINVNM);

        for (var IRHS = 1; IRHS <= NNS; IRHS++) {
          final NRHS = NSVAL[IRHS];

          // Generate NRHS random solution vectors.

          var IX = 1;
          for (var J = 1; J <= NRHS; J++) {
            zlarnv(2, ISEED, N, XACT(IX));
            IX += LDA;
          }

          for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
            // Do first for UPLO = 'U', then for UPLO = 'L'.

            final UPLO = UPLOS[IUPLO - 1];

            // Set the right hand side.

            zlaptm(UPLO, N, NRHS, ONE, D, E, XACT.asMatrix(), LDA, ZERO,
                B.asMatrix(), LDA);

            // +    TEST 2
            // Solve A*x = b and compute the residual.

            zlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);
            zpttrs(UPLO, N, NRHS, D(N + 1), E(N + 1), X.asMatrix(), LDA, INFO);

            // Check error code from ZPTTRS.

            if (INFO.value != 0) {
              alaerh(PATH, 'ZPTTRS', INFO.value, 0, UPLO, N, N, -1, -1, NRHS,
                  IMAT, NFAIL, NERRS, NOUT);
            }

            zlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
            zptt02(UPLO, N, NRHS, D, E, X.asMatrix(), LDA, WORK.asMatrix(), LDA,
                RESULT(2));

            // +    TEST 3
            // Check solution from generated exact solution.

            zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                RESULT(3));

            // +    TESTS 4, 5, and 6
            // Use iterative refinement to improve the solution.

            srnamc.SRNAMT = 'ZPTRFS';
            zptrfs(
                UPLO,
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
                RWORK,
                RWORK(NRHS + 1),
                WORK,
                RWORK(2 * NRHS + 1),
                INFO);

            // Check error code from ZPTRFS.

            if (INFO.value != 0) {
              alaerh(PATH, 'ZPTRFS', INFO.value, 0, UPLO, N, N, -1, -1, NRHS,
                  IMAT, NFAIL, NERRS, NOUT);
            }

            zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                RESULT(4));
            zptt05(N, NRHS, D, E, B.asMatrix(), LDA, X.asMatrix(), LDA,
                XACT.asMatrix(), LDA, RWORK, RWORK(NRHS + 1), RESULT(5));

            // Print information about the tests that did not pass the
            // threshold.

            for (var K = 2; K <= 6; K++) {
              if (RESULT[K] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                NOUT.println(
                    ' UPLO = \'${UPLO.a1}\', N =${N.i5}, NRHS =${NRHS.i3}, type ${IMAT.i2}, test ${K.i2}, ratio = ${RESULT[K].g12_5}');
                NFAIL++;
              }
            }
            NRUN += 5;
          }
        }
      }

      // +    TEST 7
      // Estimate the reciprocal of the condition number of the
      // matrix.

      final RCOND = Box(ZERO);
      srnamc.SRNAMT = 'ZPTCON';
      zptcon(N, D(N + 1), E(N + 1), ANORM, RCOND, RWORK, INFO);

      // Check error code from ZPTCON.

      if (INFO.value != 0) {
        alaerh(PATH, 'ZPTCON', INFO.value, 0, ' ', N, N, -1, -1, -1, IMAT,
            NFAIL, NERRS, NOUT);
      }

      RESULT[7] = dget06(RCOND.value, RCONDC);

      // Print the test ratio if greater than or equal to THRESH.

      if (RESULT[7] >= THRESH) {
        if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
        NOUT.print9999(N, IMAT, 7, RESULT[7]);
        NFAIL++;
      }
      NRUN++;
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}

extension on Nout {
  void print9999(int n, int type, int test, double ratio) {
    println(
        ' N =${n.i5}, type ${type.i2}, test ${test.i2}, ratio = ${ratio.g12_5}');
  }
}
