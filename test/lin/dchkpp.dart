import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlansp.dart';
import 'package:lapack/src/dppcon.dart';
import 'package:lapack/src/dpprfs.dart';
import 'package:lapack/src/dpptrf.dart';
import 'package:lapack/src/dpptri.dart';
import 'package:lapack/src/dpptrs.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../matgen/dlatms.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'derrpox.dart';
import 'dget04.dart';
import 'dget06.dart';
import 'dlarhs.dart';
import 'dlatb4.dart';
import 'dppt01.dart';
import 'dppt02.dart';
import 'dppt03.dart';
import 'dppt05.dart';

void dchkpp(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NNS,
  final Array<int> NSVAL_,
  final double THRESH,
  final bool TSTERR,
  final int NMAX,
  final Array<double> A_,
  final Array<double> AFAC_,
  final Array<double> AINV_,
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
  final NSVAL = NSVAL_.having();
  final A = A_.having();
  final AFAC = AFAC_.having();
  final AINV = AINV_.having();
  final B = B_.having();
  final X = X_.having();
  final XACT = XACT_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();

  const ZERO = 0.0;
  const NTYPES = 9;
  const NTESTS = 8;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'], PACKS = ['C', 'R'];
  final INFO = Box(0);
  final RCONDC = Box(0.0);

  // Initialize constants and the random number seed.

  final PATH = '${'Double precision'[0]}PP';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) derrpo(PATH, NOUT);
  infoc.INFOT = 0;

  // Do for each value of N in NVAL

  for (var IN = 1; IN <= NN; IN++) {
    final N = NVAL[IN];
    final LDA = max(N, 1);
    final XTYPE = 'N';
    final NIMAT = N <= 0 ? 1 : NTYPES;

    for (var IMAT = 1; IMAT <= NIMAT; IMAT++) {
      // Do the tests only if DOTYPE( IMAT ) is true.

      if (!DOTYPE[IMAT]) continue;

      // Skip types 3, 4, or 5 if the matrix size is too small.

      final ZEROT = IMAT >= 3 && IMAT <= 5;
      if (ZEROT && N < IMAT - 2) continue;

      // Do first for UPLO = 'U', then for UPLO = 'L'

      for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
        final UPLO = UPLOS[IUPLO - 1];
        final PACKIT = PACKS[IUPLO - 1];

        // Set up parameters with DLATB4 and generate a test matrix
        // with DLATMS.

        final (:TYPE, :KL, :KU, :ANORM, :MODE, :COND, :DIST) =
            dlatb4(PATH, IMAT, N, N);

        srnamc.SRNAMT = 'DLATMS';
        dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, COND, ANORM, KL, KU,
            PACKIT, A.asMatrix(), LDA, WORK, INFO);

        // Check error code from DLATMS.

        if (INFO.value != 0) {
          alaerh(PATH, 'DLATMS', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
              NFAIL, NERRS, NOUT);
          continue;
        }

        // For types 3-5, zero one row and column of the matrix to
        // test that INFO is returned correctly.

        final int IZERO;
        if (ZEROT) {
          if (IMAT == 3) {
            IZERO = 1;
          } else if (IMAT == 4) {
            IZERO = N;
          } else {
            IZERO = N ~/ 2 + 1;
          }

          // Set row and column IZERO of A to 0.

          if (IUPLO == 1) {
            var IOFF = (IZERO - 1) * IZERO ~/ 2;
            for (var I = 1; I <= IZERO - 1; I++) {
              A[IOFF + I] = ZERO;
            }
            IOFF += IZERO;
            for (var I = IZERO; I <= N; I++) {
              A[IOFF] = ZERO;
              IOFF += I;
            }
          } else {
            var IOFF = IZERO;
            for (var I = 1; I <= IZERO - 1; I++) {
              A[IOFF] = ZERO;
              IOFF += N - I;
            }
            IOFF -= IZERO;
            for (var I = IZERO; I <= N; I++) {
              A[IOFF + I] = ZERO;
            }
          }
        } else {
          IZERO = 0;
        }

        // Compute the L*L' or U'*U factorization of the matrix.

        final NPP = N * (N + 1) ~/ 2;
        dcopy(NPP, A, 1, AFAC, 1);
        srnamc.SRNAMT = 'DPPTRF';
        dpptrf(UPLO, N, AFAC, INFO);

        // Check error code from DPPTRF.

        if (INFO.value != IZERO) {
          alaerh(PATH, 'DPPTRF', INFO.value, IZERO, UPLO, N, N, -1, -1, -1,
              IMAT, NFAIL, NERRS, NOUT);
          continue;
        }

        // Skip the tests if INFO is not 0.

        if (INFO.value != 0) continue;

        // +    TEST 1
        // Reconstruct matrix from factors and compute residual.

        dcopy(NPP, AFAC, 1, AINV, 1);
        dppt01(UPLO, N, A, AINV, RWORK, RESULT(1));

        // +    TEST 2
        // Form the inverse and compute the residual.

        dcopy(NPP, AFAC, 1, AINV, 1);
        srnamc.SRNAMT = 'DPPTRI';
        dpptri(UPLO, N, AINV, INFO);

        // Check error code from DPPTRI.

        if (INFO.value != 0) {
          alaerh(PATH, 'DPPTRI', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
              NFAIL, NERRS, NOUT);
        }

        dppt03(
            UPLO, N, A, AINV, WORK.asMatrix(), LDA, RWORK, RCONDC, RESULT(2));

        // Print information about the tests that did not pass
        // the threshold.

        for (var K = 1; K <= 2; K++) {
          if (RESULT[K] >= THRESH) {
            if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
            NOUT.print9999(UPLO, N, IMAT, K, RESULT[K]);
            NFAIL++;
          }
        }
        NRUN += 2;

        for (var IRHS = 1; IRHS <= NNS; IRHS++) {
          final NRHS = NSVAL[IRHS];

          // +    TEST 3
          // Solve and compute residual for  A * X = B.

          srnamc.SRNAMT = 'DLARHS';
          dlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A.asMatrix(), LDA,
              XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
          dlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

          srnamc.SRNAMT = 'DPPTRS';
          dpptrs(UPLO, N, NRHS, AFAC, X.asMatrix(), LDA, INFO);

          // Check error code from DPPTRS.

          if (INFO.value != 0) {
            alaerh(PATH, 'DPPTRS', INFO.value, 0, UPLO, N, N, -1, -1, NRHS,
                IMAT, NFAIL, NERRS, NOUT);
          }

          dlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
          dppt02(UPLO, N, NRHS, A, X.asMatrix(), LDA, WORK.asMatrix(), LDA,
              RWORK, RESULT(3));

          // +    TEST 4
          // Check solution from generated exact solution.

          dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC.value,
              RESULT(4));

          // +    TESTS 5, 6, and 7
          // Use iterative refinement to improve the solution.

          srnamc.SRNAMT = 'DPPRFS';
          dpprfs(UPLO, N, NRHS, A, AFAC, B.asMatrix(), LDA, X.asMatrix(), LDA,
              RWORK, RWORK(NRHS + 1), WORK, IWORK, INFO);

          // Check error code from DPPRFS.

          if (INFO.value != 0) {
            alaerh(PATH, 'DPPRFS', INFO.value, 0, UPLO, N, N, -1, -1, NRHS,
                IMAT, NFAIL, NERRS, NOUT);
          }

          dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC.value,
              RESULT(5));
          dppt05(UPLO, N, NRHS, A, B.asMatrix(), LDA, X.asMatrix(), LDA,
              XACT.asMatrix(), LDA, RWORK, RWORK(NRHS + 1), RESULT(6));

          // Print information about the tests that did not pass
          // the threshold.

          for (var K = 3; K <= 7; K++) {
            if (RESULT[K] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.print9998(UPLO, N, NRHS, IMAT, K, RESULT[K]);
              NFAIL++;
            }
          }
          NRUN += 5;
        }

        // +    TEST 8
        // Get an estimate of RCOND = 1/COND.

        final ANORM2 = dlansp('1', UPLO, N, A, RWORK);
        final RCOND = Box(0.0);
        srnamc.SRNAMT = 'DPPCON';
        dppcon(UPLO, N, AFAC, ANORM2, RCOND, WORK, IWORK, INFO);

        // Check error code from DPPCON.

        if (INFO.value != 0) {
          alaerh(PATH, 'DPPCON', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
              NFAIL, NERRS, NOUT);
        }

        RESULT[8] = dget06(RCOND.value, RCONDC.value);

        // Print the test ratio if greater than or equal to THRESH.

        if (RESULT[8] >= THRESH) {
          if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
          NOUT.print9999(UPLO, N, IMAT, 8, RESULT[8]);
          NFAIL++;
        }
        NRUN++;
      }
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}

extension on Nout {
  void print9999(String uplo, int n, int type, int test, double ratio) {
    println(
        ' UPLO = \'${uplo.a1}\', N =${n.i5}, type ${type.i2}, test ${test.i2}, ratio =${ratio.g12_5}');
  }

  void print9998(
      String uplo, int n, int nrhs, int type, int test, double ratio) {
    println(
        ' UPLO = \'${uplo.a1}\', N =${n.i5}, NRHS=${nrhs.i3}, type ${type.i2}, test(${test.i2}) =${ratio.g12_5}');
  }
}
