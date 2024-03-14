import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zcgesv.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlaset.dart';

import '../matgen/zlatms.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'common.dart';
import 'zget08.dart';
import 'zlarhs.dart';
import 'zlatb4.dart';

void zdrvab(
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
  final Array<int> IWORK_,
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
  final IWORK = IWORK_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NTYPES = 11, NTESTS = 1;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [2006, 2007, 2008, 2009];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Zomplex precision'[0]}GE';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  infoc.INFOT = 0;

  // Do for each value of M in MVAL

  for (var IM = 1; IM <= NM; IM++) {
    final M = MVAL[IM];
    final LDA = max(1, M);

    final N = M;
    final NIMAT = M <= 0 || N <= 0 ? 1 : NTYPES;

    for (var IMAT = 1; IMAT <= NIMAT; IMAT++) {
      // Do the tests only if DOTYPE( IMAT ) is true.

      if (!DOTYPE[IMAT]) continue;

      // Skip types 5, 6, or 7 if the matrix size is too small.

      final ZEROT = IMAT >= 5 && IMAT <= 7;
      if (ZEROT && N < IMAT - 4) continue;

      // Set up parameters with ZLATB4 and generate a test matrix
      // with ZLATMS.

      final (:TYPE, :KL, :KU, :ANORM, :MODE, :CNDNUM, :DIST) =
          zlatb4(PATH, IMAT, M, N);

      srnamc.SRNAMT = 'ZLATMS';
      zlatms(M, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU,
          'No packing', A.asMatrix(), LDA, WORK, INFO);

      // Check error code from ZLATMS.

      if (INFO.value != 0) {
        alaerh(PATH, 'ZLATMS', INFO.value, 0, ' ', M, N, -1, -1, -1, IMAT,
            NFAIL, NERRS, NOUT);
        continue;
      }

      // For types 5-7, zero one or more columns of the matrix to
      // test that INFO is returned correctly.

      final int IZERO;
      if (ZEROT) {
        if (IMAT == 5) {
          IZERO = 1;
        } else if (IMAT == 6) {
          IZERO = min(M, N);
        } else {
          IZERO = min(M, N) ~/ 2 + 1;
        }
        var IOFF = (IZERO - 1) * LDA;
        if (IMAT < 7) {
          for (var I = 1; I <= M; I++) {
            A[IOFF + I] = Complex.zero;
          }
        } else {
          zlaset('Full', M, N - IZERO + 1, Complex.zero, Complex.zero,
              A(IOFF + 1).asMatrix(), LDA);
        }
      } else {
        IZERO = 0;
      }

      for (var IRHS = 1; IRHS <= NNS; IRHS++) {
        final NRHS = NSVAL[IRHS];
        final XTYPE = 'N';
        final TRANS = 'N';

        srnamc.SRNAMT = 'ZLARHS';
        zlarhs(PATH, XTYPE, ' ', TRANS, N, N, KL, KU, NRHS, A.asMatrix(), LDA,
            X.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);

        srnamc.SRNAMT = 'ZCGESV';

        zlacpy('Full', M, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);

        final ITER = Box(0);
        zcgesv(N, NRHS, A.asMatrix(), LDA, IWORK, B.asMatrix(), LDA,
            X.asMatrix(), LDA, WORK.asMatrix(), SWORK, RWORK, ITER, INFO);

        if (ITER.value < 0) {
          zlacpy('Full', M, N, AFAC.asMatrix(), LDA, A.asMatrix(), LDA);
        }

        // Check error code from ZCGESV. This should be the same as
        // the one of DGETRF.

        if (INFO.value != IZERO) {
          if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
          NERRS.value++;

          if (INFO.value != IZERO && IZERO != 0) {
            NOUT.println(
                ' *** ZCGESV returned with INFO =${INFO.value.i5} instead of ${IZERO.i5}\n ==> M =${M.i5}, type ${IMAT.i2}');
          } else {
            NOUT.println(
                ' *** Error code from ZCGESV=${INFO.value.i5} for M=${M.i5}, type ${IMAT.i2}');
          }
        }

        // Skip the remaining test if the matrix is singular.

        if (INFO.value != 0) continue;

        // Check the quality of the solution

        zlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);

        zget08(TRANS, N, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
            WORK.asMatrix(), LDA, RWORK, RESULT(1));

        // Check if the test passes the testing.
        // Print information about the tests that did not
        // pass the testing.

        // If iterative refinement has been used and claimed to
        // be successful (ITER.value>0), we want
        //   NORMI(B - A*X)/(NORMI(A)*NORMI(X)*EPS*SRQT(N)) < 1

        // If double precision has been used (ITER.value<0), we want
        //   NORMI(B - A*X)/(NORMI(A)*NORMI(X)*EPS) < THRES
        // (Cf. the linear solver testing routines)

        if ((THRESH <= 0.0e+00) ||
            ((ITER.value >= 0) &&
                (N > 0) &&
                (RESULT[1] >= sqrt(N.toDouble()))) ||
            ((ITER.value < 0) && (RESULT[1] >= THRESH))) {
          if (NFAIL == 0 && NERRS.value == 0) {
            NOUT.println('\n DGE:  General dense matrices');
            NOUT.println(' Matrix types:');
            NOUT.println(
                '    1. Diagonal${' ' * 24}7. Last n/2 columns zero\n    2. Upper triangular${' ' * 16}8. Random, CNDNUM = sqrt(0.1/EPS)\n    3. Lower triangular${' ' * 16}9. Random, CNDNUM = 0.1/EPS\n    4. Random, CNDNUM = 2${' ' * 13}10. Scaled near underflow\n    5. First column zero${' ' * 14}11. Scaled near overflow\n    6. Last column zero');
            NOUT.println(' Test ratios:');
            NOUT.println(
                '   ${1.i2}: norm_1( B - A * X )  / ( norm_1(A) * norm_1(X) * EPS * sqrt(N) ) > 1 if ITERREF\n    or norm_1( B - A * X )  / ( norm_1(A) * norm_1(X) * EPS ) > THRES if DGETRF');
            NOUT.println(' Messages:');
          }

          NOUT.println(
              ' TRANS=\'${TRANS.a1}\', N =${N.i5}, NRHS=${NRHS.i3}, type ${IMAT.i2}, test(${1.i2}) =${RESULT[1].g12_5}');
          NFAIL++;
        }
        NRUN++;
      }
    }
  }

  // Print a summary of the results.

  if (NFAIL > 0) {
    NOUT.println(
        ' ZCGESV: ${NFAIL.i6} out of ${NRUN.i6} tests failed to pass the threshold');
  } else {
    NOUT.println(
        '\n All tests for ZCGESV routines passed the threshold ( ${NRUN.i6} tests run)');
  }
  if (NERRS.value > 0) {
    NOUT.println('${' ' * 6}${NERRS.value.i6} error messages recorded');
  }
}
