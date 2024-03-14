import 'dart:math';

import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlanhb.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zpbcon.dart';
import 'package:lapack/src/zpbrfs.dart';
import 'package:lapack/src/zpbtrf.dart';
import 'package:lapack/src/zpbtrs.dart';

import '../matgen/zlatms.dart';
import 'alaerh.dart';
import 'alahd.dart';
import 'alasum.dart';
import 'common.dart';
import 'dget06.dart';
import 'xlaenv.dart';
import 'zerrpox.dart';
import 'zget04.dart';
import 'zlaipd.dart';
import 'zlarhs.dart';
import 'zlatb4.dart';
import 'zpbt01.dart';
import 'zpbt02.dart';
import 'zpbt05.dart';

void zchkpb(
  final Array<bool> DOTYPE_,
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
  final Nout NOUT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DOTYPE = DOTYPE_.having();
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

  const ONE = 1.0, ZERO = 0.0;
  const NTYPES = 8, NTESTS = 7;
  const NBW = 4;
  final ISEED = Array<int>(4), KDVAL = Array<int>(NBW);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  final PATH = '${'Zomplex precision'[0]}PB';
  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  // Test the error exits

  if (TSTERR) zerrpo(PATH, NOUT);
  infoc.INFOT = 0;
  KDVAL[1] = 0;

  // Do for each value of N in NVAL

  for (var IN = 1; IN <= NN; IN++) {
    final N = NVAL[IN];
    final LDA = max(N, 1);
    var XTYPE = 'N';

    // Set limits on the number of loop iterations.

    final NKD = max(1, min(N, 4));
    final NIMAT = N == 0 ? 1 : NTYPES;

    KDVAL[2] = N + (N + 1) ~/ 4;
    KDVAL[3] = (3 * N - 1) ~/ 4;
    KDVAL[4] = (N + 1) ~/ 4;

    for (var IKD = 1; IKD <= NKD; IKD++) {
      // Do for KD = 0, (5*N+1)/4, (3N-1)/4, and (N+1)/4. This order
      // makes it easier to skip redundant values for small values
      // of N.

      final KD = KDVAL[IKD];
      final LDAB = KD + 1;

      // Do first for UPLO = 'U', then for UPLO = 'L'

      for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
        final (UPLO, KOFF, PACKIT) =
            IUPLO == 1 ? ('U', max(1, KD + 2 - N), 'Q') : ('L', 1, 'B');

        var IZERO = 0, I1 = 0, I2 = 0;
        for (var IMAT = 1; IMAT <= NIMAT; IMAT++) {
          // Do the tests only if DOTYPE( IMAT ) is true.

          if (!DOTYPE[IMAT]) continue;

          // Skip types 2, 3, or 4 if the matrix size is too small.

          final ZEROT = IMAT >= 2 && IMAT <= 4;
          if (ZEROT && N < IMAT - 1) continue;

          if (!ZEROT || !DOTYPE[1]) {
            // Set up parameters with ZLATB4 and generate a test
            // matrix with ZLATMS.

            final (:TYPE, KL: _, KU: _, :ANORM, :MODE, :CNDNUM, :DIST) =
                zlatb4(PATH, IMAT, N, N);

            srnamc.SRNAMT = 'ZLATMS';
            zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KD, KD,
                PACKIT, A(KOFF).asMatrix(), LDAB, WORK, INFO);

            // Check error code from ZLATMS.

            if (INFO.value != 0) {
              alaerh(PATH, 'ZLATMS', INFO.value, 0, UPLO, N, N, KD, KD, -1,
                  IMAT, NFAIL, NERRS, NOUT);
              continue;
            }
          } else if (IZERO > 0) {
            // Use the same matrix for types 3 and 4 as for type
            // 2 by copying back the zeroed out column,

            var IW = 2 * LDA + 1;
            if (IUPLO == 1) {
              var IOFF = (IZERO - 1) * LDAB + KD + 1;
              zcopy(IZERO - I1, WORK(IW), 1, A(IOFF - IZERO + I1), 1);
              IW += IZERO - I1;
              zcopy(I2 - IZERO + 1, WORK(IW), 1, A(IOFF), max(LDAB - 1, 1));
            } else {
              var IOFF = (I1 - 1) * LDAB + 1;
              zcopy(IZERO - I1, WORK(IW), 1, A(IOFF + IZERO - I1),
                  max(LDAB - 1, 1));
              IOFF = (IZERO - 1) * LDAB + 1;
              IW += IZERO - I1;
              zcopy(I2 - IZERO + 1, WORK(IW), 1, A(IOFF), 1);
            }
          }

          // For types 2-4, zero one row and column of the matrix
          // to test that INFO.value is returned correctly.

          IZERO = 0;
          if (ZEROT) {
            if (IMAT == 2) {
              IZERO = 1;
            } else if (IMAT == 3) {
              IZERO = N;
            } else {
              IZERO = N ~/ 2 + 1;
            }

            // Save the zeroed out row and column in WORK(*,3)

            var IW = 2 * LDA;
            for (var I = 1; I <= min(2 * KD + 1, N); I++) {
              WORK[IW + I] = Complex.zero;
            }
            IW++;
            final I1 = max(IZERO - KD, 1);
            final I2 = min(IZERO + KD, N);

            if (IUPLO == 1) {
              var IOFF = (IZERO - 1) * LDAB + KD + 1;
              zswap(IZERO - I1, A(IOFF - IZERO + I1), 1, WORK(IW), 1);
              IW += IZERO - I1;
              zswap(I2 - IZERO + 1, A(IOFF), max(LDAB - 1, 1), WORK(IW), 1);
            } else {
              var IOFF = (I1 - 1) * LDAB + 1;
              zswap(IZERO - I1, A(IOFF + IZERO - I1), max(LDAB - 1, 1),
                  WORK(IW), 1);
              IOFF = (IZERO - 1) * LDAB + 1;
              IW += IZERO - I1;
              zswap(I2 - IZERO + 1, A(IOFF), 1, WORK(IW), 1);
            }
          }

          // Set the imaginary part of the diagonals.

          if (IUPLO == 1) {
            zlaipd(N, A(KD + 1), LDAB, 0);
          } else {
            zlaipd(N, A(1), LDAB, 0);
          }

          // Do for each value of NB in NBVAL

          for (var INB = 1; INB <= NNB; INB++) {
            final NB = NBVAL[INB];
            xlaenv(1, NB);

            // Compute the L*L' or U'*U factorization of the band
            // matrix.

            zlacpy(
                'Full', KD + 1, N, A.asMatrix(), LDAB, AFAC.asMatrix(), LDAB);
            srnamc.SRNAMT = 'ZPBTRF';
            zpbtrf(UPLO, N, KD, AFAC.asMatrix(), LDAB, INFO);

            // Check error code from ZPBTRF.

            if (INFO.value != IZERO) {
              alaerh(PATH, 'ZPBTRF', INFO.value, IZERO, UPLO, N, N, KD, KD, NB,
                  IMAT, NFAIL, NERRS, NOUT);
              continue;
            }

            // Skip the tests if INFO.value is not 0.

            if (INFO.value != 0) continue;

            // +    TEST 1
            // Reconstruct matrix from factors and compute
            // residual.

            zlacpy('Full', KD + 1, N, AFAC.asMatrix(), LDAB, AINV.asMatrix(),
                LDAB);
            zpbt01(UPLO, N, KD, A.asMatrix(), LDAB, AINV.asMatrix(), LDAB,
                RWORK, RESULT(1));

            // Print the test ratio if it is >= THRESH.

            if (RESULT[1] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
              NOUT.println(
                  ' UPLO=\'${UPLO.a1}\', N=${N.i5}, KD=${KD.i5}, NB=${NB.i4}, type ${IMAT.i2}, test ${1.i2}, ratio= ${RESULT[1].g12_5}');
              NFAIL++;
            }
            NRUN++;

            // Only do other tests if this is the first blocksize.

            if (INB > 1) continue;

            // Form the inverse of A so we can get a good estimate
            // of RCONDC = 1/(norm(A) * norm(inv(A))).

            zlaset(
                'Full', N, N, Complex.zero, Complex.one, AINV.asMatrix(), LDA);
            srnamc.SRNAMT = 'ZPBTRS';
            zpbtrs(UPLO, N, KD, N, AFAC.asMatrix(), LDAB, AINV.asMatrix(), LDA,
                INFO);

            // Compute RCONDC = 1/(norm(A) * norm(inv(A))).

            final ANORM = zlanhb('1', UPLO, N, KD, A.asMatrix(), LDAB, RWORK);
            final AINVNM = zlange('1', N, N, AINV.asMatrix(), LDA, RWORK);
            final double RCONDC;
            if (ANORM <= ZERO || AINVNM <= ZERO) {
              RCONDC = ONE;
            } else {
              RCONDC = (ONE / ANORM) / AINVNM;
            }

            for (var IRHS = 1; IRHS <= NNS; IRHS++) {
              final NRHS = NSVAL[IRHS];

              // +    TEST 2
              // Solve and compute residual for A * X = B.

              srnamc.SRNAMT = 'ZLARHS';
              zlarhs(PATH, XTYPE, UPLO, ' ', N, N, KD, KD, NRHS, A.asMatrix(),
                  LDAB, XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
              zlacpy('Full', N, NRHS, B.asMatrix(), LDA, X.asMatrix(), LDA);

              srnamc.SRNAMT = 'ZPBTRS';
              zpbtrs(UPLO, N, KD, NRHS, AFAC.asMatrix(), LDAB, X.asMatrix(),
                  LDA, INFO);

              // Check error code from ZPBTRS.

              if (INFO.value != 0) {
                alaerh(PATH, 'ZPBTRS', INFO.value, 0, UPLO, N, N, KD, KD, NRHS,
                    IMAT, NFAIL, NERRS, NOUT);
              }

              zlacpy('Full', N, NRHS, B.asMatrix(), LDA, WORK.asMatrix(), LDA);
              zpbt02(UPLO, N, KD, NRHS, A.asMatrix(), LDAB, X.asMatrix(), LDA,
                  WORK.asMatrix(), LDA, RWORK, RESULT(2));

              // +    TEST 3
              // Check solution from generated exact solution.

              zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                  RESULT(3));

              // +    TESTS 4, 5, and 6
              // Use iterative refinement to improve the solution.

              srnamc.SRNAMT = 'ZPBRFS';
              zpbrfs(
                  UPLO,
                  N,
                  KD,
                  NRHS,
                  A.asMatrix(),
                  LDAB,
                  AFAC.asMatrix(),
                  LDAB,
                  B.asMatrix(),
                  LDA,
                  X.asMatrix(),
                  LDA,
                  RWORK,
                  RWORK(NRHS + 1),
                  WORK,
                  RWORK(2 * NRHS + 1),
                  INFO);

              // Check error code from ZPBRFS.

              if (INFO.value != 0) {
                alaerh(PATH, 'ZPBRFS', INFO.value, 0, UPLO, N, N, KD, KD, NRHS,
                    IMAT, NFAIL, NERRS, NOUT);
              }

              zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA, RCONDC,
                  RESULT(4));
              zpbt05(
                  UPLO,
                  N,
                  KD,
                  NRHS,
                  A.asMatrix(),
                  LDAB,
                  B.asMatrix(),
                  LDA,
                  X.asMatrix(),
                  LDA,
                  XACT.asMatrix(),
                  LDA,
                  RWORK,
                  RWORK(NRHS + 1),
                  RESULT(5));

              // Print information about the tests that did not
              // pass the threshold.

              for (var K = 2; K <= 6; K++) {
                if (RESULT[K] >= THRESH) {
                  if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                  NOUT.println(
                      ' UPLO=\'${UPLO.a1}\', N=${N.i5}, KD=${KD.i5}, NRHS=${NRHS.i3}, type ${IMAT.i2}, test(${K.i2}) = ${RESULT[K].g12_5}');
                  NFAIL++;
                }
              }
              NRUN += 5;
            }

            // +    TEST 7
            // Get an estimate of RCOND = 1/CNDNUM.
            {
              final RCOND = Box(ZERO), ANORM = Box(ZERO);
              srnamc.SRNAMT = 'ZPBCON';
              zpbcon(UPLO, N, KD, AFAC.asMatrix(), LDAB, ANORM, RCOND, WORK,
                  RWORK, INFO);

              // Check error code from ZPBCON.

              if (INFO.value != 0) {
                alaerh(PATH, 'ZPBCON', INFO.value, 0, UPLO, N, N, KD, KD, -1,
                    IMAT, NFAIL, NERRS, NOUT);
              }

              RESULT[7] = dget06(RCOND.value, RCONDC);

              // Print the test ratio if it is >= THRESH.

              if (RESULT[7] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) alahd(NOUT, PATH);
                NOUT.println(
                    ' UPLO=\'${UPLO.a1}\', N=${N.i5}, KD=${KD.i5},${' ' * 10} type ${IMAT.i2}, test(${7.i2}) = ${RESULT[7].g12_5}');
                NFAIL++;
              }
              NRUN++;
            }
          }
        }
      }
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, NERRS.value);
}
