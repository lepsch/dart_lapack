import 'dart:math';

import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zgbcon.dart';
import 'package:lapack/src/zgbrfs.dart';
import 'package:lapack/src/zgbtrf.dart';
import 'package:lapack/src/zgbtrs.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlangb.dart';
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
import 'zgbt01.dart';
import 'zgbt02.dart';
import 'zgbt05.dart';
import 'zget04.dart';
import 'zlarhs.dart';
import 'zlatb4.dart';

void zchkgb(
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
  final Array<Complex> A_,
  final int LA,
  final Array<Complex> AFAC_,
  final int LAFAC,
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
  final B = B_.having();
  final X = X_.having();
  final XACT = XACT_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();

  const ONE = 1.0, ZERO = 0.0;
  const NTYPES = 8, NTESTS = 7;
  const NBW = 4, NTRAN = 3;
  final ISEED = Array<int>(4), KLVAL = Array<int>(NBW), KUVAL = Array<int>(NBW);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991], TRANSS = ['N', 'T', 'C'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Zomplex precision'[0]}GB';
  var NRUN = 0;
  var NFAIL = 0;
  var NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) zerrge(PATH, NOUT);
  infoc.INFOT = 0;

  // Initialize the first value for the lower and upper bandwidths.

  KLVAL[1] = 0;
  KUVAL[1] = 0;

  // Do for each value of M in MVAL

  for (var IM = 1; IM <= NM; IM++) {
    final M = MVAL[IM];

    // Set values to use for the lower bandwidth.

    KLVAL[2] = M + (M + 1) ~/ 4;

    // KLVAL( 2 ) = max( M-1, 0 )

    KLVAL[3] = (3 * M - 1) ~/ 4;
    KLVAL[4] = (M + 1) ~/ 4;

    // Do for each value of N in NVAL

    for (var IN = 1; IN <= NN; IN++) {
      final N = NVAL[IN];
      var XTYPE = 'N';

      // Set values to use for the upper bandwidth.

      KUVAL[2] = N + (N + 1) ~/ 4;

      // KUVAL( 2 ) = max( N-1, 0 )

      KUVAL[3] = (3 * N - 1) ~/ 4;
      KUVAL[4] = (N + 1) ~/ 4;

      // Set limits on the number of loop iterations.

      final NKL = N == 0 ? 2 : min(M + 1, 4);
      final NKU = M == 0 ? 2 : min(N + 1, 4);
      final NIMAT = M <= 0 || N <= 0 ? 1 : NTYPES;

      for (var IKL = 1; IKL <= NKL; IKL++) {
        // Do for KL = 0, (5*M+1)/4, (3M-1)/4, and (M+1)/4. This
        // order makes it easier to skip redundant values for small
        // values of M.

        final KL = KLVAL[IKL];
        for (var IKU = 1; IKU <= NKU; IKU++) {
          // Do for KU = 0, (5*N+1)/4, (3N-1)/4, and (N+1)/4. This
          // order makes it easier to skip redundant values for
          // small values of N.

          final KU = KUVAL[IKU];

          // Check that A and AFAC are big enough to generate this
          // matrix.

          final LDA = KL + KU + 1;
          final LDAFAC = 2 * KL + KU + 1;
          if ((LDA * N) > LA || (LDAFAC * N) > LAFAC) {
            if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
            if (N * (KL + KU + 1) > LA) {
              NOUT.println(
                  ' *** In ZCHKGB, LA=${LA.i5} is too small for M=${M.i5}, N=${N.i5}, KL=${KL.i4}, KU=${KU.i4}\n ==> Increase LA to at least ${(N * (KL + KU + 1)).i5}');
              NERRS.value += 1;
            }
            if (N * (2 * KL + KU + 1) > LAFAC) {
              NOUT.println(
                  ' *** In ZCHKGB, LAFAC=${LAFAC.i5} is too small for M=${M.i5}, N=${N.i5}, KL=${KL.i4}, KU=${KU.i4}\n ==> Increase LAFAC to at least ${(N * (2 * KL + KU + 1)).i5}');
              NERRS.value += 1;
            }
            continue;
          }

          int IZERO = 0, I1 = 0, I2 = 0, IOFF = 0;
          for (var IMAT = 1; IMAT <= NIMAT; IMAT++) {
            // Do the tests only if DOTYPE( IMAT ) is true.

            if (!DOTYPE[IMAT]) continue;

            // Skip types 2, 3, or 4 if the matrix size is too
            // small.

            final ZEROT = IMAT >= 2 && IMAT <= 4;
            if (ZEROT && N < IMAT - 1) continue;

            if (!ZEROT || !DOTYPE[1]) {
              // Set up parameters with ZLATB4 and generate a
              // test matrix with ZLATMS.

              final (:TYPE, :KL, :KU, :ANORM, :MODE, :CNDNUM, :DIST) =
                  zlatb4(PATH, IMAT, M, N);

              final KOFF = max(1, KU + 2 - N);
              for (var I = 1; I <= KOFF - 1; I++) {
                A[I] = Complex.zero;
              }
              srnamc.SRNAMT = 'ZLATMS';
              zlatms(M, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL,
                  KU, 'Z', A(KOFF).asMatrix(), LDA, WORK, INFO);

              // Check the error code from ZLATMS.

              if (INFO.value != 0) {
                alaerh(PATH, 'ZLATMS', INFO.value, 0, ' ', M, N, KL, KU, -1,
                    IMAT, NFAIL, NERRS, NOUT);
                continue;
              }
            } else if (IZERO > 0) {
              // Use the same matrix for types 3 and 4 as for
              // type 2 by copying back the zeroed out column.

              zcopy(I2 - I1 + 1, B, 1, A(IOFF + I1), 1);
            }

            // For types 2, 3, and 4, zero one or more columns of
            // the matrix to test that INFO.value is returned correctly.

            IZERO = 0;
            if (ZEROT) {
              if (IMAT == 2) {
                IZERO = 1;
              } else if (IMAT == 3) {
                IZERO = min(M, N);
              } else {
                IZERO = min(M, N) ~/ 2 + 1;
              }
              var IOFF = (IZERO - 1) * LDA;
              if (IMAT < 4) {
                // Store the column to be zeroed out in B.

                I1 = max(1, KU + 2 - IZERO);
                I2 = min(KL + KU + 1, KU + 1 + (M - IZERO));
                zcopy(I2 - I1 + 1, A(IOFF + I1), 1, B, 1);

                for (var I = I1; I <= I2; I++) {
                  A[IOFF + I] = Complex.zero;
                }
              } else {
                for (var J = IZERO; J <= N; J++) {
                  for (var I = max(1, KU + 2 - J);
                      I <= min(KL + KU + 1, KU + 1 + (M - J));
                      I++) {
                    A[IOFF + I] = Complex.zero;
                  }
                  IOFF += LDA;
                }
              }
            }

            // These lines, if used in place of the calls in the
            // loop over INB, cause the code to bomb on a Sun
            // SPARCstation.

            // ANORMO = zlangb( 'O', N, KL, KU, A, LDA, RWORK )
            // ANORMI = zlangb( 'I', N, KL, KU, A, LDA, RWORK )

            // Do for each blocksize in NBVAL

            for (var INB = 1; INB <= NNB; INB++) {
              final NB = NBVAL[INB];
              xlaenv(1, NB);

              // Compute the LU factorization of the band matrix.

              if (M > 0 && N > 0) {
                zlacpy('Full', KL + KU + 1, N, A.asMatrix(), LDA,
                    AFAC(KL + 1).asMatrix(), LDAFAC);
              }
              srnamc.SRNAMT = 'ZGBTRF';
              zgbtrf(M, N, KL, KU, AFAC.asMatrix(), LDAFAC, IWORK, INFO);

              // Check error code from ZGBTRF.

              if (INFO.value != IZERO) {
                alaerh(PATH, 'ZGBTRF', INFO.value, IZERO, ' ', M, N, KL, KU, NB,
                    IMAT, NFAIL, NERRS, NOUT);
              }

              // +    TEST 1
              // Reconstruct matrix from factors and compute
              // residual.

              zgbt01(M, N, KL, KU, A.asMatrix(), LDA, AFAC.asMatrix(), LDAFAC,
                  IWORK, WORK, RESULT(1));

              // Print information about the tests so far that
              // did not pass the threshold.

              if (RESULT[1] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                NOUT.println(
                    ' M =${M.i5}, N =${N.i5}, KL=${KL.i5}, KU=${KU.i5}, NB =${NB.i4}, type ${IMAT.i1}, test(${1.i1})=${RESULT[1].g12_5}');
                NFAIL++;
              }
              NRUN++;

              // Skip the remaining tests if this is not the
              // first block size or if M != N.

              if (INB > 1 || M != N) continue;

              final ANORMO = zlangb('O', N, KL, KU, A.asMatrix(), LDA, RWORK);
              final ANORMI = zlangb('I', N, KL, KU, A.asMatrix(), LDA, RWORK);
              final LDB = max(1, N);

              final bool TRFCON;
              final double RCONDO, RCONDI;
              if (INFO.value == 0) {
                TRFCON = false;
                // Form the inverse of A so we can get a good
                // estimate of CNDNUM = norm(A) * norm(inv(A)).

                zlaset('Full', N, N, Complex.zero, Complex.one, WORK.asMatrix(),
                    LDB);
                srnamc.SRNAMT = 'ZGBTRS';
                zgbtrs('No transpose', N, KL, KU, N, AFAC.asMatrix(), LDAFAC,
                    IWORK, WORK.asMatrix(), LDB, INFO);

                // Compute the 1-norm condition number of A.

                var AINVNM = zlange('O', N, N, WORK.asMatrix(), LDB, RWORK);
                if (ANORMO <= ZERO || AINVNM <= ZERO) {
                  RCONDO = ONE;
                } else {
                  RCONDO = (ONE / ANORMO) / AINVNM;
                }

                // Compute the infinity-norm condition number of
                // A.

                AINVNM = zlange('I', N, N, WORK.asMatrix(), LDB, RWORK);
                if (ANORMI <= ZERO || AINVNM <= ZERO) {
                  RCONDI = ONE;
                } else {
                  RCONDI = (ONE / ANORMI) / AINVNM;
                }
              } else {
                // Do only the condition estimate if INFO.value != 0.

                TRFCON = true;
                RCONDO = ZERO;
                RCONDI = ZERO;
              }

              // Skip the solve tests if the matrix is singular.

              if (!TRFCON) {
                for (var IRHS = 1; IRHS <= NNS; IRHS++) {
                  final NRHS = NSVAL[IRHS];
                  XTYPE = 'N';

                  for (var ITRAN = 1; ITRAN <= NTRAN; ITRAN++) {
                    final TRANS = TRANSS[ITRAN - 1];
                    final (RCONDC, _) =
                        ITRAN == 1 ? (RCONDO, 'O') : (RCONDI, 'I');

                    // +    TEST 2:
                    // Solve and compute residual for op(A) * X = B.

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
                        LDB,
                        B.asMatrix(),
                        LDB,
                        ISEED,
                        INFO);
                    XTYPE = 'C';
                    zlacpy(
                        'Full', N, NRHS, B.asMatrix(), LDB, X.asMatrix(), LDB);

                    srnamc.SRNAMT = 'ZGBTRS';
                    zgbtrs(TRANS, N, KL, KU, NRHS, AFAC.asMatrix(), LDAFAC,
                        IWORK, X.asMatrix(), LDB, INFO);

                    // Check error code from ZGBTRS.

                    if (INFO.value != 0) {
                      alaerh(PATH, 'ZGBTRS', INFO.value, 0, TRANS, N, N, KL, KU,
                          -1, IMAT, NFAIL, NERRS, NOUT);
                    }

                    zlacpy('Full', N, NRHS, B.asMatrix(), LDB, WORK.asMatrix(),
                        LDB);
                    zgbt02(
                        TRANS,
                        M,
                        N,
                        KL,
                        KU,
                        NRHS,
                        A.asMatrix(),
                        LDA,
                        X.asMatrix(),
                        LDB,
                        WORK.asMatrix(),
                        LDB,
                        RWORK,
                        RESULT(2));

                    // +    TEST 3:
                    // Check solution from generated exact
                    // solution.

                    zget04(N, NRHS, X.asMatrix(), LDB, XACT.asMatrix(), LDB,
                        RCONDC, RESULT(3));

                    // +    TESTS 4, 5, 6:
                    // Use iterative refinement to improve the
                    // solution.

                    srnamc.SRNAMT = 'ZGBRFS';
                    zgbrfs(
                        TRANS,
                        N,
                        KL,
                        KU,
                        NRHS,
                        A.asMatrix(),
                        LDA,
                        AFAC.asMatrix(),
                        LDAFAC,
                        IWORK,
                        B.asMatrix(),
                        LDB,
                        X.asMatrix(),
                        LDB,
                        RWORK,
                        RWORK(NRHS + 1),
                        WORK,
                        RWORK(2 * NRHS + 1),
                        INFO);

                    // Check error code from ZGBRFS.

                    if (INFO.value != 0) {
                      alaerh(PATH, 'ZGBRFS', INFO.value, 0, TRANS, N, N, KL, KU,
                          NRHS, IMAT, NFAIL, NERRS, NOUT);
                    }

                    zget04(N, NRHS, X.asMatrix(), LDB, XACT.asMatrix(), LDB,
                        RCONDC, RESULT(4));
                    zgbt05(
                        TRANS,
                        N,
                        KL,
                        KU,
                        NRHS,
                        A.asMatrix(),
                        LDA,
                        B.asMatrix(),
                        LDB,
                        X.asMatrix(),
                        LDB,
                        XACT.asMatrix(),
                        LDB,
                        RWORK,
                        RWORK(NRHS + 1),
                        RESULT(5));

                    // Print information about the tests that did
                    // not pass the threshold.

                    for (var K = 2; K <= 6; K++) {
                      if (RESULT[K] >= THRESH) {
                        if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                        NOUT.println(
                            ' TRANS=\'${TRANS.a1}\', N=${N.i5}, KL=${KL.i5}, KU=${KU.i5}, NRHS=${NRHS.i3}, type ${IMAT.i1}, test(${K.i1})=${RESULT[K].g12_5}');
                        NFAIL++;
                      }
                    }
                    NRUN += 5;
                  }
                }
              }

              // +    TEST 7:
              // Get an estimate of RCOND = 1/CNDNUM.

              for (var ITRAN = 1; ITRAN <= 2; ITRAN++) {
                final (ANORM, RCONDC, NORM) =
                    ITRAN == 1 ? (ANORMO, RCONDO, 'O') : (ANORMI, RCONDI, 'I');
                final RCOND = Box(ZERO);
                srnamc.SRNAMT = 'ZGBCON';
                zgbcon(NORM, N, KL, KU, AFAC.asMatrix(), LDAFAC, IWORK, ANORM,
                    RCOND, WORK, RWORK, INFO);

                // Check error code from ZGBCON.

                if (INFO.value != 0) {
                  alaerh(PATH, 'ZGBCON', INFO.value, 0, NORM, N, N, KL, KU, -1,
                      IMAT, NFAIL, NERRS, NOUT);
                }

                RESULT[7] = dget06(RCOND.value, RCONDC);

                // Print information about the tests that did
                // not pass the threshold.

                if (RESULT[7] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                  NOUT.println(
                      ' NORM =\'${NORM.a1}\', N=${N.i5}, KL=${KL.i5}, KU=${KU.i5},${' ' * 10} type ${IMAT.i1}, test(${7.i1})=${RESULT[7].g12_5}');
                  NFAIL++;
                }
                NRUN++;
              }
            }
          }
        }
      }
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
