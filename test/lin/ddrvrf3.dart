import 'dart:math';

import 'package:lapack/src/blas/dtrsm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgelqf.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dtfsm.dart';
import 'package:lapack/src/dtrttf.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/dgeqrf.dart';

import '../matgen/dlarnd.dart';
import 'common.dart';

void ddrvrf3(
  final Nout NOUT,
  final int NN,
  final Array<int> NVAL_,
  final double THRESH,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> ARF_,
  final Matrix<double> B1_,
  final Matrix<double> B2_,
  final Array<double> D_WORK_DLANGE_,
  final Array<double> D_WORK_DGEQRF_,
  final Array<double> TAU_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final NVAL = NVAL_.having(length: NN);
  final A = A_.having(ld: LDA);
  final ARF = ARF_.having();
  final B1 = B1_.having(ld: LDA);
  final B2 = B2_.having(ld: LDA);
  final D_WORK_DLANGE = D_WORK_DLANGE_.having();
  final D_WORK_DGEQRF = D_WORK_DGEQRF_.having();
  final TAU = TAU_.having();
  const ZERO = 0.0, ONE = 1.0;
  const NTESTS = 1;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'];
  const FORMS = ['N', 'T'];
  const SIDES = ['L', 'R'];
  const TRANSS = ['N', 'T'];
  const DIAGS = ['N', 'U'];

  // Initialize constants and the random number seed.

  var NRUN = 0;
  var NFAIL = 0;
  final INFO = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }
  final EPS = dlamch('Precision');

  for (var IIM = 1; IIM <= NN; IIM++) {
    final M = NVAL[IIM];

    for (var IIN = 1; IIN <= NN; IIN++) {
      final N = NVAL[IIN];

      for (var IFORM = 1; IFORM <= 2; IFORM++) {
        final CFORM = FORMS[IFORM - 1];

        for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
          final UPLO = UPLOS[IUPLO - 1];

          for (var ISIDE = 1; ISIDE <= 2; ISIDE++) {
            final SIDE = SIDES[ISIDE - 1];

            for (var ITRANS = 1; ITRANS <= 2; ITRANS++) {
              final TRANS = TRANSS[ITRANS - 1];

              for (var IDIAG = 1; IDIAG <= 2; IDIAG++) {
                final DIAG = DIAGS[IDIAG - 1];

                for (var IALPHA = 1; IALPHA <= 3; IALPHA++) {
                  final double ALPHA;
                  if (IALPHA == 1) {
                    ALPHA = ZERO;
                  } else if (IALPHA == 2) {
                    ALPHA = ONE;
                  } else {
                    ALPHA = dlarnd(2, ISEED);
                  }

                  // All the parameters are set:
                  //    CFORM, SIDE, UPLO, TRANS, DIAG, M, N,
                  //    and ALPHA
                  // READY TO TEST!

                  NRUN++;

                  final int NA;
                  if (ISIDE == 1) {
                    // The case ISIDE == 1 is when SIDE == 'L'
                    // -> A is M-by-M ( B is M-by-N )

                    NA = M;
                  } else {
                    // The case ISIDE == 2 is when SIDE == 'R'
                    // -> A is N-by-N ( B is M-by-N )

                    NA = N;
                  }

                  // Generate A our NA--by--NA triangular
                  // matrix.
                  // Our test is based on forward error so we
                  // do want A to be well conditioned! To get
                  // a well-conditioned triangular matrix, we
                  // take the R factor of the QR/LQ factorization
                  // of a random matrix.

                  for (var J = 1; J <= NA; J++) {
                    for (var I = 1; I <= NA; I++) {
                      A[I][J] = dlarnd(2, ISEED);
                    }
                  }

                  if (IUPLO == 1) {
                    // The case IUPLO == 1 is when SIDE == 'U'
                    // -> QR factorization.

                    srnamc.SRNAMT = 'DGEQRF';
                    dgeqrf(NA, NA, A, LDA, TAU, D_WORK_DGEQRF, LDA, INFO);

                    // Forcing main diagonal of test matrix to
                    // be unit makes it ill-conditioned for
                    // some test cases

                    if (lsame(DIAG, 'U')) {
                      for (var J = 1; J <= NA; J++) {
                        for (var I = 1; I <= J; I++) {
                          A[I][J] /= 2.0 * A[J][J];
                        }
                      }
                    }
                  } else {
                    // The case IUPLO == 2 is when SIDE == 'L'
                    // -> QL factorization.

                    srnamc.SRNAMT = 'DGELQF';
                    dgelqf(NA, NA, A, LDA, TAU, D_WORK_DGEQRF, LDA, INFO);

                    // Forcing main diagonal of test matrix to
                    // be unit makes it ill-conditioned for
                    // some test cases

                    if (lsame(DIAG, 'U')) {
                      for (var I = 1; I <= NA; I++) {
                        for (var J = 1; J <= I; J++) {
                          A[I][J] /= 2.0 * A[I][I];
                        }
                      }
                    }
                  }

                  // Store a copy of A in RFP format (in ARF).

                  srnamc.SRNAMT = 'DTRTTF';
                  dtrttf(CFORM, UPLO, NA, A, LDA, ARF, INFO);

                  // Generate B1 our M--by--N right-hand side
                  // and store a copy in B2.

                  for (var J = 1; J <= N; J++) {
                    for (var I = 1; I <= M; I++) {
                      B1[I][J] = dlarnd(2, ISEED);
                      B2[I][J] = B1[I][J];
                    }
                  }

                  // Solve op( A ) X = B or X op( A ) = B
                  // with DTRSM

                  srnamc.SRNAMT = 'DTRSM';
                  dtrsm(SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, LDA, B1, LDA);

                  // Solve op( A ) X = B or X op( A ) = B
                  // with DTFSM

                  srnamc.SRNAMT = 'DTFSM';
                  dtfsm(CFORM, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, ARF, B2,
                      LDA);

                  // Check that the result agrees.

                  for (var J = 1; J <= N; J++) {
                    for (var I = 1; I <= M; I++) {
                      B1[I][J] = B2[I][J] - B1[I][J];
                    }
                  }

                  RESULT[1] = dlange('I', M, N, B1, LDA, D_WORK_DLANGE);

                  RESULT[1] /= sqrt(EPS) / max(max(M, N), 1);

                  if (RESULT[1] >= THRESH) {
                    if (NFAIL == 0) {
                      NOUT.println();
                      NOUT.println(
                          '  *** Error(s) or Failure(s) while testing DTFSM ***');
                    }
                    NOUT.println(
                        '      Failure in DTFSM, CFORM=\'${CFORM.a1}\', SIDE=\'${SIDE.a1}\', UPLO=\'${UPLO.a1}\', TRANS=\'${TRANS.a1}\', DIAG=\'${DIAG.a1}\', M=${M.i3}, N =${N.i3}, test=${RESULT[1].g12_5}');
                    NFAIL++;
                  }
                }
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
        ' All tests for DTFSM auxiliary routine passed the threshold ( ${NRUN.i5} tests run)');
  } else {
    NOUT.println(
        ' DTFSM  auxiliary routine: ${NFAIL.i5} out of ${NRUN.i5} tests failed to pass the threshold');
  }
}
