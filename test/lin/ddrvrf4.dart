import 'dart:math';

import 'package:lapack/src/blas/dsyrk.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dsfrk.dart';
import 'package:lapack/src/dtfttr.dart';
import 'package:lapack/src/dtrttf.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../matgen/dlarnd.dart';
import 'common.dart';

void ddrvrf4(
  final Nout NOUT,
  final int NN,
  final Array<int> NVAL_,
  final double THRESH,
  final Matrix<double> C1_,
  final Matrix<double> C2_,
  final int LDC,
  final Array<double> CRF_,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> D_WORK_DLANGE_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final NVAL = NVAL_.having(length: NN);
  final C1 = C1_.having(ld: LDC);
  final C2 = C2_.having(ld: LDC);
  final CRF = CRF_.having();
  final A = A_.having(ld: LDA);
  final D_WORK_DLANGE = D_WORK_DLANGE_.having();
  const ZERO = 0.0, ONE = 1.0;
  const NTESTS = 1;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'];
  const FORMS = ['N', 'T'];
  const TRANSS = ['N', 'T'];

  // Initialize constants and the random number seed.

  var NRUN = 0;
  var NFAIL = 0;
  final INFO = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }
  final EPS = dlamch('Precision');

  for (var IIN = 1; IIN <= NN; IIN++) {
    final N = NVAL[IIN];

    for (var IIK = 1; IIK <= NN; IIK++) {
      final K = NVAL[IIN];

      for (var IFORM = 1; IFORM <= 2; IFORM++) {
        final CFORM = FORMS[IFORM - 1];

        for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
          final UPLO = UPLOS[IUPLO - 1];

          for (var ITRANS = 1; ITRANS <= 2; ITRANS++) {
            final TRANS = TRANSS[ITRANS - 1];

            for (var IALPHA = 1; IALPHA <= 4; IALPHA++) {
              final double ALPHA, BETA;
              if (IALPHA == 1) {
                ALPHA = ZERO;
                BETA = ZERO;
              } else if (IALPHA == 2) {
                ALPHA = ONE;
                BETA = ZERO;
              } else if (IALPHA == 3) {
                ALPHA = ZERO;
                BETA = ONE;
              } else {
                ALPHA = dlarnd(2, ISEED);
                BETA = dlarnd(2, ISEED);
              }

              // All the parameters are set:
              //    CFORM, UPLO, TRANS, M, N,
              //    ALPHA, and BETA
              // READY TO TEST!

              NRUN++;

              final double NORMA;
              if (ITRANS == 1) {
                // In this case we are NOTRANS, so A is N-by-K

                for (var J = 1; J <= K; J++) {
                  for (var I = 1; I <= N; I++) {
                    A[I][J] = dlarnd(2, ISEED);
                  }
                }

                NORMA = dlange('I', N, K, A, LDA, D_WORK_DLANGE);
              } else {
                // In this case we are TRANS, so A is K-by-N

                for (var J = 1; J <= N; J++) {
                  for (var I = 1; I <= K; I++) {
                    A[I][J] = dlarnd(2, ISEED);
                  }
                }

                NORMA = dlange('I', K, N, A, LDA, D_WORK_DLANGE);
              }

              // Generate C1 our N--by--N symmetric matrix.
              // Make sure C2 has the same upper/lower part,
              // (the one that we do not touch), so
              // copy the initial C1 in C2 in it.

              for (var J = 1; J <= N; J++) {
                for (var I = 1; I <= N; I++) {
                  C1[I][J] = dlarnd(2, ISEED);
                  C2[I][J] = C1[I][J];
                }
              }

              // (See comment later on for why we use DLANGE and
              // not DLANSY for C1.)

              dlange('I', N, N, C1, LDC, D_WORK_DLANGE);

              srnamc.SRNAMT = 'DTRTTF';
              dtrttf(CFORM, UPLO, N, C1, LDC, CRF, INFO);

              // call dsyrk the BLAS routine -> gives C1

              srnamc.SRNAMT = 'DSYRK ';
              dsyrk(UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C1, LDC);

              // call dsfrk the RFP routine -> gives CRF

              srnamc.SRNAMT = 'DSFRK ';
              dsfrk(CFORM, UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, CRF);

              // convert CRF in full format -> gives C2

              srnamc.SRNAMT = 'DTFTTR';
              dtfttr(CFORM, UPLO, N, CRF, C2, LDC, INFO);

              // compare C1 and C2

              for (var J = 1; J <= N; J++) {
                for (var I = 1; I <= N; I++) {
                  C1[I][J] -= C2[I][J];
                }
              }

              // Yes, C1 is symmetric so we could call DLANSY,
              // but we want to check the upper part that is
              // supposed to be unchanged and the diagonal that
              // is supposed to be real -> DLANGE

              RESULT[1] = dlange('I', N, N, C1, LDC, D_WORK_DLANGE);
              RESULT[1] = RESULT[1] /
                  max(ALPHA.abs() * NORMA + BETA.abs(), ONE) /
                  max(N, 1) /
                  EPS;

              if (RESULT[1] >= THRESH) {
                if (NFAIL == 0) {
                  NOUT.println();
                  NOUT.println(
                      '  *** Error(s) or Failure(s) while testing DSFRK ***');
                }
                NOUT.println(
                    '      Failure in DSFRK, CFORM=\'${CFORM.a1}\', UPLO=\'${UPLO.a1}\', TRANS=\'${TRANS.a1}\', N=${N.i3}, K =${K.i3}, test=${RESULT[1].g12_5}');
                NFAIL++;
              }
            }
          }
        }
      }
    }
  }

  // Print a summary of the results.

  if (NFAIL == 0) {
    NOUT.println(
        ' All tests for DSFRK auxiliary routine passed the threshold ( ${NRUN.i5} tests run)');
  } else {
    NOUT.println(
        ' DSFRK auxiliary routine: ${NFAIL.i5} out of ${NRUN.i5} tests failed to pass the threshold');
  }
}
