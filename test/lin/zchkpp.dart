import 'dart:math';

import 'package:lapack/lapack.dart';

import '../matgen/zlatms.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'dget06.dart';
import 'zerrpo.dart';
import 'zget04.dart';
import 'zlaipd.dart';
import 'zlarhs.dart';
import 'zlatb4.dart';
import 'zppt01.dart';
import 'zppt02.dart';
import 'zppt03.dart';
import 'zppt05.dart';

void zchkpp(
  final Array<bool> DOTYPE_,
  final int NN,
  final Array<int> NVAL_,
  final int NNS,
  final Array<int> NSVAL_,
  final double THRESH,
  final bool TSTERR,
  final int NMAX,
  final Array<Complex> A_,
  final Array<Complex> AFAC_,
  final Array<Complex> AINV_,
  final Array<Complex> B_,
  final Array<Complex> X_,
  final Array<Complex> XACT_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
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

  const ZERO = 0.0;
  const NTYPES = 9, NTESTS = 8;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'], PACKS = ['C', 'R'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Zomplex precision'[0]}PP';
  var NRUN = 0;
  var NFAIL = 0;
  var NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) zerrpo(PATH, NOUT);
  infoc.INFOT = 0;

  // Do for each value of N in NVAL

  for (var IN = 1; IN <= NN; IN++) {
    final N = NVAL[IN];
    final LDA = max(N, 1);
    var XTYPE = 'N';
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

        // Set up parameters with ZLATB4 and generate a test matrix
        // with ZLATMS.

        final (:TYPE, :KL, :KU, :ANORM, :MODE, :CNDNUM, :DIST) =
            zlatb4(PATH, IMAT, N, N);

        srnamc.SRNAMT = 'ZLATMS';
        zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
            PACKIT, A.asMatrix(), LDA, WORK, INFO);

        // Check error code from ZLATMS.

        if (INFO.value != 0) {
          alaerh(PATH, 'ZLATMS', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
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
              A[IOFF + I] = Complex.zero;
            }
            IOFF += IZERO;
            for (var I = IZERO; I <= N; I++) {
              A[IOFF] = Complex.zero;
              IOFF += I;
            }
          } else {
            var IOFF = IZERO;
            for (var I = 1; I <= IZERO - 1; I++) {
              A[IOFF] = Complex.zero;
              IOFF += N - I;
            }
            IOFF -= IZERO;
            for (var I = IZERO; I <= N; I++) {
              A[IOFF + I] = Complex.zero;
            }
          }
        } else {
          IZERO = 0;
        }

        // Set the imaginary part of the diagonals.

        if (IUPLO == 1) {
          zlaipd(N, A, 2, 1);
        } else {
          zlaipd(N, A, N, -1);
        }

        // Compute the L*L' or U'*U factorization of the matrix.

        final NPP = N * (N + 1) ~/ 2;
        zcopy(NPP, A, 1, AFAC, 1);
        srnamc.SRNAMT = 'ZPPTRF';
        zpptrf(UPLO, N, AFAC, INFO);

        // Check error code from ZPPTRF.

        if (INFO.value != IZERO) {
          alaerh(PATH, 'ZPPTRF', INFO.value, IZERO, UPLO, N, N, -1, -1, -1,
              IMAT, NFAIL, NERRS, NOUT);
          continue;
        }

        // Skip the tests if INFO is not 0.

        if (INFO.value != 0) continue;

        // +    TEST 1
        // Reconstruct matrix from factors and compute residual.

        zcopy(NPP, AFAC, 1, AINV, 1);
        zppt01(UPLO, N, A, AINV, RWORK, RESULT(1));

        // +    TEST 2
        // Form the inverse and compute the residual.

        zcopy(NPP, AFAC, 1, AINV, 1);
        srnamc.SRNAMT = 'ZPPTRI';
        zpptri(UPLO, N, AINV, INFO);

        // Check error code from ZPPTRI.

        if (INFO.value != 0) {
          alaerh(PATH, 'ZPPTRI', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
              NFAIL, NERRS, NOUT);
        }

        final RCONDC = Box(ZERO);
        zppt03(
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

          srnamc.SRNAMT = 'ZLARHS';
          zlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A.asMatrix(), LDA,
              XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
          zlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

          srnamc.SRNAMT = 'ZPPTRS';
          zpptrs(UPLO, N, NRHS, AFAC, X.asMatrix(), LDA, INFO);

          // Check error code from ZPPTRS.

          if (INFO.value != 0) {
            alaerh(PATH, 'ZPPTRS', INFO.value, 0, UPLO, N, N, -1, -1, NRHS,
                IMAT, NFAIL, NERRS, NOUT);
          }

          zlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
          zppt02(UPLO, N, NRHS, A, X.asMatrix(), LDA, WORK.asMatrix(), LDA,
              RWORK, RESULT(3));

          // +    TEST 4
          // Check solution from generated exact solution.

          zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC.value,
              RESULT(4));

          // +    TESTS 5, 6, and 7
          // Use iterative refinement to improve the solution.

          srnamc.SRNAMT = 'ZPPRFS';
          zpprfs(UPLO, N, NRHS, A, AFAC, B.asMatrix(), LDA, X.asMatrix(), LDA,
              RWORK, RWORK(NRHS + 1), WORK, RWORK(2 * NRHS + 1), INFO);

          // Check error code from ZPPRFS.

          if (INFO.value != 0) {
            alaerh(PATH, 'ZPPRFS', INFO.value, 0, UPLO, N, N, -1, -1, NRHS,
                IMAT, NFAIL, NERRS, NOUT);
          }

          zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC.value,
              RESULT(5));
          zppt05(UPLO, N, NRHS, A, B.asMatrix(), LDA, X.asMatrix(), LDA,
              XACT.asMatrix(), LDA, RWORK, RWORK(NRHS + 1), RESULT(6));

          // Print information about the tests that did not pass
          // the threshold.

          for (var K = 3; K <= 7; K++) {
            if (RESULT[K] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.println(
                  ' UPLO = \'${UPLO.a1}\', N =${N.i5}, NRHS=${NRHS.i3}, type ${IMAT.i2}, test(${K.i2}) =${RESULT[K].g12_5}');
              NFAIL++;
            }
          }
          NRUN += 5;
        }

        // +    TEST 8
        // Get an estimate of RCOND = 1/CNDNUM.

        final ANORM2 = zlanhp('1', UPLO, N, A, RWORK);
        final RCOND = Box(ZERO);
        srnamc.SRNAMT = 'ZPPCON';
        zppcon(UPLO, N, AFAC, ANORM2, RCOND, WORK, RWORK, INFO);

        // Check error code from ZPPCON.

        if (INFO.value != 0) {
          alaerh(PATH, 'ZPPCON', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
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
}
