import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zcposv.dart';
import 'package:lapack/src/zlacpy.dart';

import '../matgen/zlatms.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'common.dart';
import 'zlaipd.dart';
import 'zlarhs.dart';
import 'zlatb4.dart';
import 'zpot06.dart';

void zdrvac(
  final Array<bool> DOTYPE_,
  final int NM,
  final Array<int> MVAL_,
  final int NNS,
  final Array<int> NSVAL_,
  final double THRESH,
  final int NMAX,
  final Array<Complex> A_,
  final Array<Complex> AFAC_,
  final Array<Complex> B_,
  final Array<Complex> X_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Array<Complex> SWORK_,
  final Nout NOUT,
) {
  final DOTYPE = DOTYPE_.having();
  final MVAL = MVAL_.having();
  final NSVAL = NSVAL_.having();
  final A = A_.having();
  final AFAC = AFAC_.having();
  final B = B_.having();
  final X = X_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final SWORK = SWORK_.having();
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NTYPES = 9, NTESTS = 1;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Zomplex precision'[0]}PO';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  infoc.INFOT = 0;

  // Do for each value of N in MVAL

  for (var IM = 1; IM <= NM; IM++) {
    final N = MVAL[IM];
    final LDA = max(N, 1);
    final NIMAT = N <= 0 ? 1 : NTYPES;

    imatLoop:
    for (var IMAT = 1; IMAT <= NIMAT; IMAT++) {
      // Do the tests only if DOTYPE( IMAT ) is true.

      if (!DOTYPE[IMAT]) continue;

      // Skip types 3, 4, or 5 if the matrix size is too small.

      final ZEROT = IMAT >= 3 && IMAT <= 5;
      if (ZEROT && N < IMAT - 2) continue;

      // Do first for UPLO = 'U', then for UPLO = 'L'

      for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
        final UPLO = UPLOS[IUPLO - 1];

        // Set up parameters with ZLATB4 and generate a test matrix
        // with ZLATMS.

        final (:TYPE, :KL, :KU, :ANORM, :MODE, :CNDNUM, :DIST) =
            zlatb4(PATH, IMAT, N, N);

        srnamc.SRNAMT = 'ZLATMS';
        zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
            UPLO, A.asMatrix(), LDA, WORK, INFO);

        // Check error code from ZLATMS.

        if (INFO.value != 0) {
          alaerh(PATH, 'ZLATMS', INFO.value, 0, UPLO, N, N, -1, -1, -1, IMAT,
              NFAIL, NERRS, NOUT);
        } else {
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
            var IOFF = (IZERO - 1) * LDA;

            // Set row and column IZERO of A to 0.

            if (IUPLO == 1) {
              for (var I = 1; I <= IZERO - 1; I++) {
                A[IOFF + I] = Complex.zero;
              }
              IOFF += IZERO;
              for (var I = IZERO; I <= N; I++) {
                A[IOFF] = Complex.zero;
                IOFF += LDA;
              }
            } else {
              IOFF = IZERO;
              for (var I = 1; I <= IZERO - 1; I++) {
                A[IOFF] = Complex.zero;
                IOFF += LDA;
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

          zlaipd(N, A, LDA + 1, 0);

          for (var IRHS = 1; IRHS <= NNS; IRHS++) {
            final NRHS = NSVAL[IRHS];
            final XTYPE = 'N';

            // Form an exact solution and set the right hand side.

            srnamc.SRNAMT = 'ZLARHS';
            zlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A.asMatrix(),
                LDA, X.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);

            // Compute the L*L' or U'*U factorization of the
            // matrix and solve the system.

            srnamc.SRNAMT = 'ZCPOSV';

            zlacpy('All', N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);

            final ITER = Box(0);
            zcposv(UPLO, N, NRHS, AFAC.asMatrix(), LDA, B.asMatrix(), LDA,
                X.asMatrix(), LDA, WORK.asMatrix(), SWORK, RWORK, ITER, INFO);

            if (ITER.value < 0) {
              zlacpy('All', N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);
            }

            // Check error code from ZCPOSV .

            if (INFO.value != IZERO) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NERRS.value++;

              if (INFO.value != IZERO && IZERO != 0) {
                NOUT.println(
                    ' *** ZCPOSV returned with INFO =${INFO.value.i5} instead of ${IZERO.i5}\n ==> N =${N.i5}, type ${IMAT.i2}');
              } else {
                NOUT.println(
                    ' *** Error code from ZCPOSV=${INFO.value.i5} for M=${N.i5}, type ${IMAT.i2}');
              }
            }

            // Skip the remaining test if the matrix is singular.

            if (INFO.value != 0) continue imatLoop;

            // Check the quality of the solution

            zlacpy('All', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);

            zpot06(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                WORK.asMatrix(), LDA, RWORK, RESULT(1));

            // Check if the test passes the testing.
            // Print information about the tests that did not
            // pass the testing.

            // If iterative refinement has been used and claimed to
            // be successful (ITER>0), we want
            // NORM1(B - A*X)/(NORM1(A)*NORM1(X)*EPS*SRQT(N)) < 1

            // If double precision has been used (ITER<0), we want
            // NORM1(B - A*X)/(NORM1(A)*NORM1(X)*EPS) < THRES
            // (Cf. the linear solver testing routines)

            if ((THRESH <= 0.0e+00) ||
                ((ITER.value >= 0) && (N > 0) && (RESULT[1] >= sqrt(N))) ||
                ((ITER.value < 0) && (RESULT[1] >= THRESH))) {
              if (NFAIL == 0 && NERRS.value == 0) {
                NOUT.println('\n ZPO:  positive definite dense matrices');
                NOUT.println(' Matrix types:');
                NOUT.println(
                    '    1. Diagonal${' ' * 24}7. Last n/2 columns zero\n    2. Upper triangular${' ' * 16}8. Random, CNDNUM = sqrt(0.1/EPS)\n    3. Lower triangular${' ' * 16}9. Random, CNDNUM = 0.1/EPS\n    4. Random, CNDNUM = 2${' ' * 13}10. Scaled near underflow\n    5. First column zero${' ' * 14}11. Scaled near overflow\n    6. Last column zero');
                NOUT.println(' Test ratios:');
                NOUT.println(
                    '   ${1.i2}: norm_1( B - A * X )  / ( norm_1(A) * norm_1(X) * EPS * sqrt(N) ) > 1 if ITERREF\n    or norm_1( B - A * X )  / ( norm_1(A) * norm_1(X) * EPS ) > THRES if ZPOTRF');
                NOUT.println(' Messages:');
              }

              NOUT.println(
                  ' UPLO=\'${UPLO.a1}\', N =${N.i5}, NRHS=${NRHS.i3}, type ${IMAT.i2}, test(${1.i2}) =${RESULT[1].g12_5}');

              NFAIL++;
            }

            NRUN++;
          }
        }
      }
    }

    // Print a summary of the results.

    if (NFAIL > 0) {
      NOUT.println(
          ' ZCPOSV: ${NFAIL.i6} out of ${NRUN.i6} tests failed to pass the threshold');
    } else {
      NOUT.println(
          '\n All tests for ZCPOSV routines passed the threshold ( ${NRUN.i6} tests run)');
    }
    if (NERRS.value > 0) {
      NOUT.println('${' ' * 6}${NERRS.value.i6} error messages recorded');
    }
  }
}
