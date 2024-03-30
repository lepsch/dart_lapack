import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zgecon.dart';
import 'package:lapack/src/zgerfs.dart';
import 'package:lapack/src/zgetrf.dart';
import 'package:lapack/src/zgetri.dart';
import 'package:lapack/src/zgetrs.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlaset.dart';

import '../matgen/zlatms.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'dget06.dart';
import 'xlaenv.dart';
import 'zerrgex.dart';
import 'zget01.dart';
import 'zget02.dart';
import 'zget03.dart';
import 'zget04.dart';
import 'zget07.dart';
import 'zlarhs.dart';
import 'zlatb4.dart';

void zchkge(
  final Array<bool> DOTYPE_,
  final int NM,
  final Array<int> MVAL_,
  final int NN,
  final Array<int> NVAL_,
  final int NNB,
  final Array<int> NBVAL_,
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
  final Array<int> IWORK_,
  final Nout NOUT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DOTYPE = DOTYPE_.having();
  final MVAL = MVAL_.having();
  final NVAL = NVAL_.having();
  final NBVAL = NBVAL_.having();
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

  const ONE = 1.0, ZERO = 0.0;
  const NTYPES = 11;
  const NTESTS = 8;
  const NTRAN = 3;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991], TRANSS = ['N', 'T', 'C'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Zomplex precision'[0]}GE';
  var NRUN = 0;
  var NFAIL = 0;
  var NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  xlaenv(1, 1);
  if (TSTERR) zerrge(PATH, NOUT);
  infoc.INFOT = 0;
  xlaenv(2, 2);

  // Do for each value of M in MVAL

  for (var IM = 1; IM <= NM; IM++) {
    final M = MVAL[IM];
    final LDA = max(1, M);

    // Do for each value of N in NVAL

    for (var IN = 1; IN <= NN; IN++) {
      final N = NVAL[IN];
      var XTYPE = 'N';
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
        // test that INFO.value is returned correctly.

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

        // These lines, if used in place of the calls in the DO 60
        // loop, cause the code to bomb on a Sun SPARCstation.

        // ANORMO = zlange( 'O', M, N, A, LDA, RWORK )
        // ANORMI = zlange( 'I', M, N, A, LDA, RWORK )

        // Do for each blocksize in NBVAL

        for (var INB = 1; INB <= NNB; INB++) {
          final NB = NBVAL[INB];
          xlaenv(1, NB);

          // Compute the LU factorization of the matrix.

          zlacpy('Full', M, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);
          srnamc.SRNAMT = 'ZGETRF';
          zgetrf(M, N, AFAC.asMatrix(), LDA, IWORK, INFO);

          // Check error code from ZGETRF.

          if (INFO.value != IZERO) {
            alaerh(PATH, 'ZGETRF', INFO.value, IZERO, ' ', M, N, -1, -1, NB,
                IMAT, NFAIL, NERRS, NOUT);
          }

          // +    TEST 1
          // Reconstruct matrix from factors and compute residual.

          zlacpy('Full', M, N, AFAC.asMatrix(), LDA, AINV.asMatrix(), LDA);
          zget01(M, N, A.asMatrix(), LDA, AINV.asMatrix(), LDA, IWORK, RWORK,
              RESULT(1));

          // +    TEST 2
          // Form the inverse if the factorization was successful
          // and compute the residual.

          final int NT;
          final bool TRFCON;
          final double ANORMI, ANORMO, RCONDI;
          final RCONDO = Box(0.0);
          if (M == N && INFO.value == 0) {
            TRFCON = false;
            zlacpy('Full', N, N, AFAC.asMatrix(), LDA, AINV.asMatrix(), LDA);
            srnamc.SRNAMT = 'ZGETRI';
            final NRHS = NSVAL[1];
            final LWORK = NMAX * max(3, NRHS).toInt();
            zgetri(N, AINV.asMatrix(), LDA, IWORK, WORK, LWORK, INFO);

            // Check error code from ZGETRI.

            if (INFO.value != 0) {
              alaerh(PATH, 'ZGETRI', INFO.value, 0, ' ', N, N, -1, -1, NB, IMAT,
                  NFAIL, NERRS, NOUT);
            }

            // Compute the residual for the matrix times its
            // inverse.  Also compute the 1-norm condition number
            // of A.

            zget03(N, A.asMatrix(), LDA, AINV.asMatrix(), LDA, WORK.asMatrix(),
                LDA, RWORK, RCONDO, RESULT(2));
            ANORMO = zlange('O', M, N, A.asMatrix(), LDA, RWORK);

            // Compute the infinity-norm condition number of A.

            ANORMI = zlange('I', M, N, A.asMatrix(), LDA, RWORK);
            final AINVNM = zlange('I', N, N, AINV.asMatrix(), LDA, RWORK);
            if (ANORMI <= ZERO || AINVNM <= ZERO) {
              RCONDI = ONE;
            } else {
              RCONDI = (ONE / ANORMI) / AINVNM;
            }
            NT = 2;
          } else {
            // Do only the condition estimate if INFO.value > 0.

            TRFCON = true;
            ANORMO = zlange('O', M, N, A.asMatrix(), LDA, RWORK);
            ANORMI = zlange('I', M, N, A.asMatrix(), LDA, RWORK);
            RCONDO.value = ZERO;
            RCONDI = ZERO;
            NT = 1;
          }

          // Print information about the tests so far that did not
          // pass the threshold.

          for (var K = 1; K <= NT; K++) {
            if (RESULT[K] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.println(
                  ' M = ${M.i5}, N =${N.i5}, NB =${NB.i4}, type ${IMAT.i2}, test(${K.i2}) =${RESULT[K].g12_5}');
              NFAIL++;
            }
          }
          NRUN += NT;

          // Skip the remaining tests if this is not the first
          // block size or if M != N.  Skip the solve tests if
          // the matrix is singular.

          if (INB > 1 || M != N) continue;
          if (!TRFCON) {
            for (var IRHS = 1; IRHS <= NNS; IRHS++) {
              final NRHS = NSVAL[IRHS];
              XTYPE = 'N';

              for (var ITRAN = 1; ITRAN <= NTRAN; ITRAN++) {
                final TRANS = TRANSS[ITRAN - 1];
                final RCONDC = ITRAN == 1 ? RCONDO.value : RCONDI;

                // +    TEST 3
                // Solve and compute residual for A * X = B.

                srnamc.SRNAMT = 'ZLARHS';
                zlarhs(
                    PATH,
                    XTYPE,
                    ' ',
                    TRANS,
                    N,
                    N,
                    KL,
                    KU,
                    NRHS,
                    A.asMatrix(),
                    LDA,
                    XACT.asMatrix(),
                    LDA,
                    B.asMatrix(),
                    LDA,
                    ISEED,
                    INFO);
                XTYPE = 'C';

                zlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);
                srnamc.SRNAMT = 'ZGETRS';
                zgetrs(TRANS, N, NRHS, AFAC.asMatrix(), LDA, IWORK,
                    X.asMatrix(), LDA, INFO);

                // Check error code from ZGETRS.

                if (INFO.value != 0) {
                  alaerh(PATH, 'ZGETRS', INFO.value, 0, TRANS, N, N, -1, -1,
                      NRHS, IMAT, NFAIL, NERRS, NOUT);
                }

                zlacpy(
                    'Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
                zget02(TRANS, N, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                    WORK.asMatrix(), LDA, RWORK, RESULT(3));

                // +    TEST 4
                // Check solution from generated exact solution.

                zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                    RESULT(4));

                // +    TESTS 5, 6, and 7
                // Use iterative refinement to improve the
                // solution.

                srnamc.SRNAMT = 'ZGERFS';
                zgerfs(
                    TRANS,
                    N,
                    NRHS,
                    A.asMatrix(),
                    LDA,
                    AFAC.asMatrix(),
                    LDA,
                    IWORK,
                    B.asMatrix(),
                    LDA,
                    X.asMatrix(),
                    LDA,
                    RWORK,
                    RWORK(NRHS + 1),
                    WORK,
                    RWORK(2 * NRHS + 1),
                    INFO);

                // Check error code from ZGERFS.

                if (INFO.value != 0) {
                  alaerh(PATH, 'ZGERFS', INFO.value, 0, TRANS, N, N, -1, -1,
                      NRHS, IMAT, NFAIL, NERRS, NOUT);
                }

                zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                    RESULT(5));
                zget07(
                    TRANS,
                    N,
                    NRHS,
                    A.asMatrix(),
                    LDA,
                    B.asMatrix(),
                    LDA,
                    X.asMatrix(),
                    LDA,
                    XACT.asMatrix(),
                    LDA,
                    RWORK,
                    true,
                    RWORK(NRHS + 1),
                    RESULT(6));

                // Print information about the tests that did not
                // pass the threshold.

                for (var K = 3; K <= 7; K++) {
                  if (RESULT[K] >= THRESH) {
                    if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                    NOUT.println(
                        ' TRANS=\'${TRANS.a1}\', N =${N.i5}, NRHS=${NRHS.i3}, type ${IMAT.i2}, test(${K.i2}) =${RESULT[K].g12_5}');
                    NFAIL++;
                  }
                }
                NRUN += 5;
              }
            }
          }

          // +    TEST 8
          // Get an estimate of RCOND = 1/CNDNUM.

          for (var ITRAN = 1; ITRAN <= 2; ITRAN++) {
            final (ANORM, RCONDC, NORM) = ITRAN == 1
                ? (ANORMO, RCONDO.value, 'O')
                : (ANORMI, RCONDI, 'I');

            final RCOND = Box(0.0);
            srnamc.SRNAMT = 'ZGECON';
            zgecon(
                NORM, N, AFAC.asMatrix(), LDA, ANORM, RCOND, WORK, RWORK, INFO);

            // Check error code from ZGECON.

            if (INFO.value != 0) {
              alaerh(PATH, 'ZGECON', INFO.value, 0, NORM, N, N, -1, -1, -1,
                  IMAT, NFAIL, NERRS, NOUT);
            }

            RESULT[8] = dget06(RCOND.value, RCONDC);

            // Print information about the tests that did not pass
            // the threshold.

            if (RESULT[8] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.println(
                  ' NORM =\'${NORM.a1}\', N =${N.i5},${' ' * 10} type ${IMAT.i2}, test(${8.i2}) =${RESULT[8].g12_5}');
              NFAIL++;
            }
            NRUN++;
          }
        }
      }
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
