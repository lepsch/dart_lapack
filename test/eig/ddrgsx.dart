import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgesvd.dart';
import 'package:lapack/src/dggesx.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/range.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:test/test.dart';

import '../matgen/dlakf2.dart';
import '../matgen/dlatm5.dart';
import '../test_driver.dart';
import 'alasvm.dart';
import 'common.dart';
import 'dget51.dart';
import 'dget53.dart';
import 'dlctsx.dart';

Future<void> ddrgsx(
  final int NSIZE,
  final int NCMAX,
  final double THRESH,
  final Nin NIN,
  final Nout NOUT,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final Matrix<double> AI_,
  final Matrix<double> BI_,
  final Matrix<double> Z_,
  final Matrix<double> Q_,
  final Array<double> ALPHAR_,
  final Array<double> ALPHAI_,
  final Array<double> BETA_,
  final Matrix<double> C_,
  final int LDC,
  final Array<double> S_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Array<bool> BWORK_,
  final Box<int> INFO,
  final TestDriver test,
  final String group,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDA);
  final AI = AI_.having(ld: LDA);
  final BI = BI_.having(ld: LDA);
  final Z = Z_.having(ld: LDA);
  final Q = Q_.having(ld: LDA);
  final ALPHAR = ALPHAR_.having();
  final ALPHAI = ALPHAI_.having();
  final BETA = BETA_.having();
  final C = C_.having(ld: LDC);
  final S = S_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  final BWORK = BWORK_.having();
  const ZERO = 0.0, ONE = 1.0, TEN = 1.0e+1;
  final RESULT = Array<double>(10);

  // Check for errors
  var MAXWRK = 0;
  {
    if (NSIZE < 0) {
      INFO.value = -1;
    } else if (THRESH < ZERO) {
      INFO.value = -2;
      // } else if ( NIN == 0 ) {
      //    INFO.value = -3;
      // } else if ( NOUT == 0 ) {
      //    INFO.value = -4;
    } else if (LDA < 1 || LDA < NSIZE) {
      INFO.value = -6;
    } else if (LDC < 1 || LDC < NSIZE * NSIZE ~/ 2) {
      INFO.value = -17;
    } else if (LIWORK < NSIZE + 6) {
      INFO.value = -21;
    }

    // Compute workspace
    // (Note: Comments in the code beginning "Workspace:" describe the
    // minimal amount of workspace needed at that point in the code,
    // as well as the preferred amount for good performance.
    // NB refers to the optimal block size for the immediately
    // following subroutine, as returned by ILAENV.)

    var MINWRK = 1;
    if (INFO.value == 0 && LWORK >= 1) {
      MINWRK = max(10 * (NSIZE + 1), 5 * NSIZE * NSIZE ~/ 2);

      // workspace for sggesx

      MAXWRK = 9 * (NSIZE + 1) +
          NSIZE * ilaenv(1, 'DGEQRF', ' ', NSIZE, 1, NSIZE, 0);
      MAXWRK = max(
          MAXWRK,
          9 * (NSIZE + 1) +
              NSIZE * ilaenv(1, 'DORGQR', ' ', NSIZE, 1, NSIZE, -1));

      // workspace for dgesvd

      final BDSPAC = 5 * NSIZE * NSIZE ~/ 2;
      MAXWRK = max(
          MAXWRK,
          3 * NSIZE * NSIZE ~/ 2 +
              NSIZE *
                  NSIZE *
                  ilaenv(
                    1,
                    'DGEBRD',
                    ' ',
                    NSIZE * NSIZE ~/ 2,
                    NSIZE * NSIZE ~/ 2,
                    -1,
                    -1,
                  ));
      MAXWRK = max(MAXWRK, BDSPAC);

      MAXWRK = max(MAXWRK, MINWRK);

      WORK[1] = MAXWRK.toDouble();
    }

    if (LWORK < MINWRK) INFO.value = -19;

    if (INFO.value != 0) {
      xerbla('DDRGSX', -INFO.value);
      return;
    }
  }

  // Important constants
  final ULP = dlamch('P');
  final ULPINV = ONE / ULP;
  final SMLNUM = dlamch('S') / ULP;
  final THRSH2 = TEN * THRESH;

  var NTESTT = 0;
  var NERRS = 0;
  final PARAMS = claenv.IPARMS.copy();

  // Go to the tests for read-in matrix pairs

  if (NSIZE != 0) {
    test.group(group, () {
      test.setUp(() {
        claenv.IPARMS.assign(PARAMS);
      });

      // Test the built-in matrix pairs.
      // Loop over different functions (IFUNC) of DGGESX, types (PRTYPE)
      // of test matrices, different size (M+N)

      var WEIGHT = sqrt(ULP);
      for (final IFUNC in 0.through(3)) {
        for (final PRTYPE in 1.through(5)) {
          for (final M in 1.through(NSIZE - 1)) {
            for (final N in 1.through(NSIZE - M)) {
              WEIGHT = ONE / WEIGHT;
              final ctx = (WEIGHT: WEIGHT);

              test('DDRGSX (FUNC=$IFUNC, PRTYPE=$PRTYPE, M=$M, N=$N)', () {
                final (:WEIGHT) = ctx;
                mn.M = M;
                mn.N = N;
                final MPLUSN = mn.MPLUSN = M + N;

                // Generate test matrices

                mn.FS = true;
                mn.K = 0;

                dlaset('Full', MPLUSN, MPLUSN, ZERO, ZERO, AI, LDA);
                dlaset('Full', MPLUSN, MPLUSN, ZERO, ZERO, BI, LDA);

                final QBA = Box(3), QBB = Box(4);
                dlatm5(
                    PRTYPE,
                    M,
                    N,
                    AI,
                    LDA,
                    AI(M + 1, M + 1),
                    LDA,
                    AI(1, M + 1),
                    LDA,
                    BI,
                    LDA,
                    BI(M + 1, M + 1),
                    LDA,
                    BI(1, M + 1),
                    LDA,
                    Q,
                    LDA,
                    Z,
                    LDA,
                    WEIGHT,
                    QBA,
                    QBB);

                // Compute the Schur factorization and swapping the
                // m-by-m (1,1)-blocks with n-by-n (2,2)-blocks.
                // Swapping is accomplished via the function dlctsx
                // which is supplied below.

                final SENSE = switch (IFUNC) {
                  0 => 'N',
                  1 => 'E',
                  2 => 'V',
                  3 => 'B',
                  _ => throw UnimplementedError(),
                };

                dlacpy('Full', MPLUSN, MPLUSN, AI, LDA, A, LDA);
                dlacpy('Full', MPLUSN, MPLUSN, BI, LDA, B, LDA);

                final PL = Array<double>(2), DIFEST = Array<double>(2);
                final MM = Box(0), LINFO = Box(0);
                dggesx(
                    'V',
                    'V',
                    'S',
                    dlctsx,
                    SENSE,
                    MPLUSN,
                    AI,
                    LDA,
                    BI,
                    LDA,
                    MM,
                    ALPHAR,
                    ALPHAI,
                    BETA,
                    Q,
                    LDA,
                    Z,
                    LDA,
                    PL,
                    DIFEST,
                    WORK,
                    LWORK,
                    IWORK,
                    LIWORK,
                    BWORK,
                    LINFO);

                if (LINFO.value != 0 && LINFO.value != MPLUSN + 2) {
                  RESULT[1] = ULPINV;
                  NOUT.println(
                      ' DDRGSX: DGGESX returned INFO=${LINFO.value.i6}.\n${' ' * 9}N=${MPLUSN.i6}, JTYPE=${PRTYPE.i6})');
                  INFO.value = LINFO.value;
                  return;
                }

                // Compute the norm(A, B)

                dlacpy('Full', MPLUSN, MPLUSN, AI, LDA, WORK.asMatrix(MPLUSN),
                    MPLUSN);
                dlacpy('Full', MPLUSN, MPLUSN, BI, LDA,
                    WORK(MPLUSN * MPLUSN + 1).asMatrix(MPLUSN), MPLUSN);
                final ABNRM = dlange('Fro', MPLUSN, 2 * MPLUSN,
                    WORK.asMatrix(MPLUSN), MPLUSN, WORK);

                // Do tests (1) to (4)

                dget51(1, MPLUSN, A, LDA, AI, LDA, Q, LDA, Z, LDA, WORK,
                    RESULT(1));
                dget51(1, MPLUSN, B, LDA, BI, LDA, Q, LDA, Z, LDA, WORK,
                    RESULT(2));
                dget51(3, MPLUSN, B, LDA, BI, LDA, Q, LDA, Q, LDA, WORK,
                    RESULT(3));
                dget51(3, MPLUSN, B, LDA, BI, LDA, Z, LDA, Z, LDA, WORK,
                    RESULT(4));
                var NTEST = 4;

                // Do tests (5) and (6): check Schur form of A and
                // compare eigenvalues with diagonals.

                var TEMP1 = ZERO;
                RESULT[5] = ZERO;
                RESULT[6] = ZERO;

                for (var J = 1; J <= MPLUSN; J++) {
                  var ILABAD = false;
                  final TEMP2 = Box(ZERO);
                  if (ALPHAI[J] == ZERO) {
                    TEMP2.value = ((ALPHAR[J] - AI[J][J]).abs() /
                                max(
                                  max(SMLNUM, ALPHAR[J].abs()),
                                  AI[J][J].abs(),
                                ) +
                            (BETA[J] - BI[J][J]).abs() /
                                max(
                                  SMLNUM,
                                  max(BETA[J].abs(), BI[J][J].abs()),
                                )) /
                        ULP;
                    if (J < MPLUSN) {
                      if (AI[J + 1][J] != ZERO) {
                        ILABAD = true;
                        RESULT[5] = ULPINV;
                      }
                    }
                    if (J > 1) {
                      if (AI[J][J - 1] != ZERO) {
                        ILABAD = true;
                        RESULT[5] = ULPINV;
                      }
                    }
                  } else {
                    final int I1;
                    if (ALPHAI[J] > ZERO) {
                      I1 = J;
                    } else {
                      I1 = J - 1;
                    }
                    if (I1 <= 0 || I1 >= MPLUSN) {
                      ILABAD = true;
                    } else if (I1 < MPLUSN - 1) {
                      if (AI[I1 + 2][I1 + 1] != ZERO) {
                        ILABAD = true;
                        RESULT[5] = ULPINV;
                      }
                    } else if (I1 > 1) {
                      if (AI[I1][I1 - 1] != ZERO) {
                        ILABAD = true;
                        RESULT[5] = ULPINV;
                      }
                    }
                    if (!ILABAD) {
                      final IINFO = Box(0);
                      dget53(AI(I1, I1), LDA, BI(I1, I1), LDA, BETA[J],
                          ALPHAR[J], ALPHAI[J], TEMP2, IINFO);
                      if (IINFO.value >= 3) {
                        print9997(NOUT, IINFO.value, J, MPLUSN, PRTYPE);
                        INFO.value = (IINFO.value).abs();
                      }
                    } else {
                      TEMP2.value = ULPINV;
                    }
                  }
                  TEMP1 = max(TEMP1, TEMP2.value);
                  if (ILABAD) {
                    print9996(NOUT, J, MPLUSN, PRTYPE);
                  }
                }
                RESULT[6] = TEMP1;
                NTEST += 2;

                // Test (7) (if sorting worked)

                RESULT[7] = ZERO;
                if (LINFO.value == MPLUSN + 3) {
                  RESULT[7] = ULPINV;
                } else if (MM.value != N) {
                  RESULT[7] = ULPINV;
                }
                NTEST++;

                // Test (8): compare the estimated value DIF and its
                // value. first, compute the exact DIF.

                RESULT[8] = ZERO;
                final MN2 = MM.value * (MPLUSN - MM.value) * 2;
                double DIFTRU = ZERO;
                if (IFUNC >= 2 && MN2 <= NCMAX * NCMAX) {
                  // Note: for either following two causes, there are
                  // almost same number of test cases fail the test.

                  dlakf2(
                      MM.value,
                      MPLUSN - MM.value,
                      AI,
                      LDA,
                      AI(MM.value + 1, MM.value + 1),
                      BI,
                      BI(MM.value + 1, MM.value + 1),
                      C,
                      LDC);

                  dgesvd('N', 'N', MN2, MN2, C, LDC, S, WORK.asMatrix(1), 1,
                      WORK(2).asMatrix(1), 1, WORK(3), LWORK - 2, INFO);
                  DIFTRU = S[MN2];

                  if (DIFEST[2] == ZERO) {
                    if (DIFTRU > ABNRM * ULP) RESULT[8] = ULPINV;
                  } else if (DIFTRU == ZERO) {
                    if (DIFEST[2] > ABNRM * ULP) RESULT[8] = ULPINV;
                  } else if ((DIFTRU > THRSH2 * DIFEST[2]) ||
                      (DIFTRU * THRSH2 < DIFEST[2])) {
                    RESULT[8] = max(DIFTRU / DIFEST[2], DIFEST[2] / DIFTRU);
                  }
                  NTEST++;
                }

                // Test (9)

                RESULT[9] = ZERO;
                if (LINFO.value == (MPLUSN + 2)) {
                  if (DIFTRU > ABNRM * ULP) RESULT[9] = ULPINV;
                  if ((IFUNC > 1) && (DIFEST[2] != ZERO)) RESULT[9] = ULPINV;
                  if ((IFUNC == 1) && (PL[1] != ZERO)) RESULT[9] = ULPINV;
                  NTEST++;
                }

                NTESTT += NTEST;

                // Print out tests which fail.

                for (var J = 1; J <= 9; J++) {
                  final reason = RESULT[J] < 10000.0
                      ? ' Matrix order=${MPLUSN.i2}, type=${PRTYPE.i2}, a=${WEIGHT.d10_3}, order(A_11)=${M.i2}, result ${J.i2} is ${RESULT[J].f8_2}'
                      : ' Matrix order=${MPLUSN.i2}, type=${PRTYPE.i2}, a=${WEIGHT.d10_3}, order(A_11)=${M.i2}, result ${J.i2} is ${RESULT[J].d10_3}';
                  test.expect(RESULT[J], lessThan(THRESH), reason: reason);
                  if (RESULT[J] >= THRESH) {
                    // If this is the first test to fail,
                    // print a header to the data file.

                    if (NERRS == 0) {
                      print9995(NOUT);

                      // Matrix types

                      NOUT.println(
                          ' Matrix types: \n  1:  A is a block diagonal matrix of Jordan blocks and B is the identity \n      matrix, \n  2:  A and B are upper triangular matrices, \n  3:  A and B are as type 2, but each second diagonal block in A_11 and \n      each third diagonal block in A_22 are 2x2 blocks,\n  4:  A and B are block diagonal matrices, \n  5:  (A,B) has potentially close or common eigenvalues.\n');

                      // Tests performed

                      print9992(NOUT);
                    }
                    NERRS++;
                    NOUT.println(reason);
                  }
                }
              });
            }
          }
        }
      }
    });
  } else {
    // Read in data from file to check accuracy of condition estimation
    // Read input data until N=0

    var NPTKNT = 0;

    while (true) {
      final double PLTRU, DIFTRU;
      final int MPLUSN, N;
      try {
        MPLUSN = mn.MPLUSN = await NIN.readInt();
        if (MPLUSN == 0) break;
        N = N = await NIN.readInt();
        await NIN.readMatrix(AI, MPLUSN, MPLUSN);
        await NIN.readMatrix(BI, MPLUSN, MPLUSN);
        (PLTRU, DIFTRU) = await NIN.readDouble2();
      } on EOF catch (_) {
        break;
      }

      final ctx = (NPTKNT: ++NPTKNT, AI: AI.copy(), BI: BI.copy());

      test.group(group, () {
        test.setUp(() {
          claenv.IPARMS.assign(PARAMS);
        });

        test('DDRGSX - precomputed (M=$MPLUSN, N=$N)', () {
          final (:NPTKNT, :AI, :BI) = ctx;
          mn.N = N;
          mn.FS = true;
          mn.K = 0;
          mn.MPLUSN = MPLUSN;
          mn.M = MPLUSN - N;

          dlacpy('Full', MPLUSN, MPLUSN, AI, LDA, A, LDA);
          dlacpy('Full', MPLUSN, MPLUSN, BI, LDA, B, LDA);

          // Compute the Schur factorization while swapping the
          // m-by-m (1,1)-blocks with n-by-n (2,2)-blocks.

          final PL = Array<double>(2), DIFEST = Array<double>(2);
          final MM = Box(0), LINFO = Box(0);
          dggesx(
              'V',
              'V',
              'S',
              dlctsx,
              'B',
              MPLUSN,
              AI,
              LDA,
              BI,
              LDA,
              MM,
              ALPHAR,
              ALPHAI,
              BETA,
              Q,
              LDA,
              Z,
              LDA,
              PL,
              DIFEST,
              WORK,
              LWORK,
              IWORK,
              LIWORK,
              BWORK,
              LINFO);

          if (LINFO.value != 0 && LINFO.value != MPLUSN + 2) {
            RESULT[1] = ULPINV;
            NOUT.println(
                ' DDRGSX: DGGESX returned INFO=${LINFO.value.i6}.\n${' ' * 9}N=${MPLUSN.i6}, Input Example #${NPTKNT.i2})');
            return;
          }

          // Compute the norm(A, B)
          // (should this be norm of (A,B) or (AI,BI)?)

          dlacpy(
              'Full', MPLUSN, MPLUSN, AI, LDA, WORK.asMatrix(MPLUSN), MPLUSN);
          dlacpy('Full', MPLUSN, MPLUSN, BI, LDA,
              WORK(MPLUSN * MPLUSN + 1).asMatrix(MPLUSN), MPLUSN);
          final ABNRM = dlange(
              'Fro', MPLUSN, 2 * MPLUSN, WORK.asMatrix(MPLUSN), MPLUSN, WORK);

          // Do tests (1) to (4)

          dget51(1, MPLUSN, A, LDA, AI, LDA, Q, LDA, Z, LDA, WORK, RESULT(1));
          dget51(1, MPLUSN, B, LDA, BI, LDA, Q, LDA, Z, LDA, WORK, RESULT(2));
          dget51(3, MPLUSN, B, LDA, BI, LDA, Q, LDA, Q, LDA, WORK, RESULT(3));
          dget51(3, MPLUSN, B, LDA, BI, LDA, Z, LDA, Z, LDA, WORK, RESULT(4));

          // Do tests (5) and (6): check Schur form of A and compare
          // eigenvalues with diagonals.

          var NTEST = 6;
          var TEMP1 = ZERO;
          RESULT[5] = ZERO;
          RESULT[6] = ZERO;

          for (var J = 1; J <= MPLUSN; J++) {
            var ILABAD = false;
            final TEMP2 = Box(ZERO);
            if (ALPHAI[J] == ZERO) {
              TEMP2.value = ((ALPHAR[J] - AI[J][J]).abs() /
                          max(max(SMLNUM, ALPHAR[J].abs()), AI[J][J].abs()) +
                      (BETA[J] - BI[J][J]).abs() /
                          max(SMLNUM, max(BETA[J].abs(), BI[J][J].abs()))) /
                  ULP;
              if (J < MPLUSN) {
                if (AI[J + 1][J] != ZERO) {
                  ILABAD = true;
                  RESULT[5] = ULPINV;
                }
              }
              if (J > 1) {
                if (AI[J][J - 1] != ZERO) {
                  ILABAD = true;
                  RESULT[5] = ULPINV;
                }
              }
            } else {
              final int I1;
              if (ALPHAI[J] > ZERO) {
                I1 = J;
              } else {
                I1 = J - 1;
              }
              if (I1 <= 0 || I1 >= MPLUSN) {
                ILABAD = true;
              } else if (I1 < MPLUSN - 1) {
                if (AI[I1 + 2][I1 + 1] != ZERO) {
                  ILABAD = true;
                  RESULT[5] = ULPINV;
                }
              } else if (I1 > 1) {
                if (AI[I1][I1 - 1] != ZERO) {
                  ILABAD = true;
                  RESULT[5] = ULPINV;
                }
              }
              if (!ILABAD) {
                final IINFO = Box(0);
                dget53(AI(I1, I1), LDA, BI(I1, I1), LDA, BETA[J], ALPHAR[J],
                    ALPHAI[J], TEMP2, IINFO);
                if (IINFO.value >= 3) {
                  print9997(NOUT, IINFO.value, J, MPLUSN, NPTKNT);
                  INFO.value = (IINFO.value).abs();
                }
              } else {
                TEMP2.value = ULPINV;
              }
            }
            TEMP1 = max(TEMP1, TEMP2.value);
            if (ILABAD) {
              print9996(NOUT, J, MPLUSN, NPTKNT);
            }
          }
          RESULT[6] = TEMP1;

          // Test (7) (if sorting worked)  <--------- need to be checked.

          NTEST = 7;
          RESULT[7] = ZERO;
          if (LINFO.value == MPLUSN + 3) RESULT[7] = ULPINV;

          // Test (8): compare the estimated value of DIF and its true value.

          NTEST = 8;
          RESULT[8] = ZERO;
          if (DIFEST[2] == ZERO) {
            if (DIFTRU > ABNRM * ULP) RESULT[8] = ULPINV;
          } else if (DIFTRU == ZERO) {
            if (DIFEST[2] > ABNRM * ULP) RESULT[8] = ULPINV;
          } else if ((DIFTRU > THRSH2 * DIFEST[2]) ||
              (DIFTRU * THRSH2 < DIFEST[2])) {
            RESULT[8] = max(DIFTRU / DIFEST[2], DIFEST[2] / DIFTRU);
          }

          // Test (9)

          NTEST = 9;
          RESULT[9] = ZERO;
          if (LINFO.value == (MPLUSN + 2)) {
            if (DIFTRU > ABNRM * ULP) RESULT[9] = ULPINV;
            // if ((IFUNC > 1) && (DIFEST[2] != ZERO)) RESULT[9] = ULPINV;
            // if ((IFUNC == 1) && (PL[1] != ZERO)) RESULT[9] = ULPINV;
          }

          // Test (10): compare the estimated value of PL and it true value.

          NTEST = 10;
          RESULT[10] = ZERO;
          if (PL[1] == ZERO) {
            if (PLTRU > ABNRM * ULP) RESULT[10] = ULPINV;
          } else if (PLTRU == ZERO) {
            if (PL[1] > ABNRM * ULP) RESULT[10] = ULPINV;
          } else if ((PLTRU > THRESH * PL[1]) || (PLTRU * THRESH < PL[1])) {
            RESULT[10] = ULPINV;
          }

          NTESTT += NTEST;

          // Print out tests which fail.

          for (var J = 1; J <= NTEST; J++) {
            final reason = RESULT[J] < 10000.0
                ? ' Input example #${NPTKNT.i2}, matrix order=${MPLUSN.i4}, result ${J.i2} is${RESULT[J].f8_2}'
                : ' Input example #${NPTKNT.i2}, matrix order=${MPLUSN.i4}, result ${J.i2} is${(RESULT[J] * 10).d10_3}';
            test.expect(RESULT[J], lessThan(THRESH), reason: reason);
            if (RESULT[J] >= THRESH) {
              // If this is the first test to fail,
              // print a header to the data file.

              if (NERRS == 0) {
                print9995(NOUT);

                // Matrix types

                NOUT.println('Input Example');

                // Tests performed

                print9992(NOUT);
              }
              NERRS++;
              NOUT.println(reason);
            }
          }
        });
      });
    }
  }

  // Summary
  alasvm('DGX', NOUT, NERRS, NTESTT, 0);

  WORK[1] = MAXWRK.toDouble();
}

void print9997(
  final Nout nout,
  final int info,
  final int eigenvalue,
  final int n,
  final int jtype,
) {
  nout.println(
      ' DDRGSX: DGET53 returned INFO=${info.i1} for eigenvalue ${eigenvalue.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6})');
}

void print9996(
  final Nout nout,
  final int eigenvalue,
  final int n,
  final int jtype,
) {
  nout.println(
      ' DDRGSX: S not in Schur form at eigenvalue ${eigenvalue.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6})');
}

void print9995(final Nout nout) {
  nout.println('\n DGX -- Real Expert Generalized Schur form problem driver');
}

void print9992(final Nout nout) {
  nout.println(
      '\n Tests performed:  (S is Schur, T is triangular, Q and Z are orthogonal,\n${' ' * 19} a is alpha, b is beta, and \' means transpose.)\n  1 = | A - Q S Z\' | / ( |A| n ulp )      2 = | B - Q T Z\' | / ( |B| n ulp )\n  3 = | I - QQ\' | / ( n ulp )             4 = | I - ZZ\' | / ( n ulp )\n  5 = 1/ULP  if A is not in Schur form S\n  6 = difference between (alpha,beta) and diagonals of (S,T)\n  7 = 1/ULP  if SDIM is not the correct number of selected eigenvalues\n  8 = 1/ULP  if DIFEST/DIFTRU > 10*THRESH or DIFTRU/DIFEST > 10*THRESH\n  9 = 1/ULP  if DIFEST <> 0 or DIFTRU > ULP*norm(A,B) when reordering fails\n 10 = 1/ULP  if PLEST/PLTRU > THRESH or PLTRU/PLEST > THRESH\n    ( Test 10 is only for input examples )\n');
}
