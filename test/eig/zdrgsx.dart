import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgesvd.dart';
import 'package:lapack/src/zggesx.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlaset.dart';

import '../matgen/zlakf2.dart';
import '../matgen/zlatm5.dart';
import 'alasvm.dart';
import 'common.dart';
import 'zget51.dart';
import 'zlctsx.dart';

Future<void> zdrgsx(
  final int NSIZE,
  final int NCMAX,
  final double THRESH,
  final Nin NIN,
  final Nout NOUT,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final Matrix<Complex> AI_,
  final Matrix<Complex> BI_,
  final Matrix<Complex> Z_,
  final Matrix<Complex> Q_,
  final Array<Complex> ALPHA_,
  final Array<Complex> BETA_,
  final Matrix<Complex> C_,
  final int LDC,
  final Array<double> S_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final int LIWORK,
  final Array<bool> BWORK_,
  final Box<int> INFO,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final C = C_.having(ld: LDC);
  final B = B_.having(ld: LDA);
  final AI = AI_.having(ld: LDA);
  final BI = BI_.having(ld: LDA);
  final Z = Z_.having(ld: LDA);
  final Q = Q_.having(ld: LDA);
  final ALPHA = ALPHA_.having();
  final BETA = BETA_.having();
  final S = S_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  final BWORK = BWORK_.having();
  const ZERO = 0.0, ONE = 1.0, TEN = 1.0e+1;
  bool ILABAD;
  String SENSE = '';
  int BDSPAC,
      IFUNC,
      J,
      MAXWRK = 0,
      MINWRK,
      MN2,
      NERRS,
      NPTKNT,
      NTEST,
      NTESTT,
      PRTYPE;
  double ABNRM,
      DIFTRU = 0,
      PLTRU = 0,
      SMLNUM,
      TEMP1,
      TEMP2,
      THRSH2,
      ULP,
      ULPINV,
      WEIGHT;
  final DIFEST = Array<double>(2),
      PL = Array<double>(2),
      RESULT = Array<double>(10);
  final LINFO = Box(0), MM = Box(0), QBA = Box(0), QBB = Box(0);

  double ABS1(Complex X) => X.toDouble().abs() + X.imaginary.abs();

  // Check for errors

  INFO.value = 0;
  if (NSIZE < 0) {
    INFO.value = -1;
  } else if (THRESH < ZERO) {
    INFO.value = -2;
    // } else if ( NIN <= 0 ) {
    //    INFO.value = -3;
    // } else if ( NOUT <= 0 ) {
    //    INFO.value = -4;
  } else if (LDA < 1 || LDA < NSIZE) {
    INFO.value = -6;
  } else if (LDC < 1 || LDC < NSIZE * NSIZE ~/ 2) {
    INFO.value = -15;
  } else if (LIWORK < NSIZE + 2) {
    INFO.value = -21;
  }

  // Compute workspace
  //  (Note: Comments in the code beginning "Workspace:" describe the
  //   minimal amount of workspace needed at that point in the code,
  //   as well as the preferred amount for good performance.
  //   NB refers to the optimal block size for the immediately
  //   following subroutine, as returned by ILAENV.)

  MINWRK = 1;
  if (INFO.value == 0 && LWORK >= 1) {
    MINWRK = 3 * NSIZE * NSIZE ~/ 2;

    // workspace for cggesx

    MAXWRK = NSIZE * (1 + ilaenv(1, 'ZGEQRF', ' ', NSIZE, 1, NSIZE, 0));
    MAXWRK = max(
        MAXWRK, NSIZE * (1 + ilaenv(1, 'ZUNGQR', ' ', NSIZE, 1, NSIZE, -1)));

    // workspace for zgesvd

    BDSPAC = 3 * NSIZE * NSIZE ~/ 2;
    MAXWRK = max(
        MAXWRK,
        NSIZE *
            NSIZE *
            (1 +
                ilaenv(1, 'ZGEBRD', ' ', NSIZE * NSIZE ~/ 2, NSIZE * NSIZE ~/ 2,
                    -1, -1)));
    MAXWRK = max(MAXWRK, BDSPAC);

    MAXWRK = max(MAXWRK, MINWRK);

    WORK[1] = MAXWRK.toComplex();
  }

  if (LWORK < MINWRK) INFO.value = -18;

  if (INFO.value != 0) {
    xerbla('ZDRGSX', -INFO.value);
    return;
  }

  // Important constants

  ULP = dlamch('P');
  ULPINV = ONE / ULP;
  SMLNUM = dlamch('S') / ULP;
  THRSH2 = TEN * THRESH;
  NTESTT = 0;
  NERRS = 0;

  // Go to the tests for read-in matrix pairs

  IFUNC = 0;
  if (NSIZE != 0) {
    // Test the built-in matrix pairs.
    // Loop over different functions (IFUNC) of ZGGESX, types (PRTYPE)
    // of test matrices, different size (mn.M+mn.N)

    PRTYPE = 0;
    QBA.value = 3;
    QBB.value = 4;
    WEIGHT = sqrt(ULP);

    for (IFUNC = 0; IFUNC <= 3; IFUNC++) {
      // 60
      for (PRTYPE = 1; PRTYPE <= 5; PRTYPE++) {
        // 50
        for (mn.M = 1; mn.M <= NSIZE - 1; mn.M++) {
          // 40
          for (mn.N = 1; mn.N <= NSIZE - mn.M; mn.N++) {
            // 30

            WEIGHT = ONE / WEIGHT;
            mn.MPLUSN = mn.M + mn.N;

            // Generate test matrices

            mn.FS = true;
            mn.K = 0;

            zlaset('Full', mn.MPLUSN, mn.MPLUSN, Complex.zero, Complex.zero, AI,
                LDA);
            zlaset('Full', mn.MPLUSN, mn.MPLUSN, Complex.zero, Complex.zero, BI,
                LDA);

            zlatm5(
                PRTYPE,
                mn.M,
                mn.N,
                AI,
                LDA,
                AI(mn.M + 1, mn.M + 1),
                LDA,
                AI(1, mn.M + 1),
                LDA,
                BI,
                LDA,
                BI(mn.M + 1, mn.M + 1),
                LDA,
                BI(1, mn.M + 1),
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
            // Swapping is accomplished via the function zlctsx
            // which is supplied below.

            if (IFUNC == 0) {
              SENSE = 'mn.N';
            } else if (IFUNC == 1) {
              SENSE = 'E';
            } else if (IFUNC == 2) {
              SENSE = 'V';
            } else if (IFUNC == 3) {
              SENSE = 'B';
            }

            zlacpy('Full', mn.MPLUSN, mn.MPLUSN, AI, LDA, A, LDA);
            zlacpy('Full', mn.MPLUSN, mn.MPLUSN, BI, LDA, B, LDA);

            zggesx(
                'V',
                'V',
                'S',
                zlctsx,
                SENSE,
                mn.MPLUSN,
                AI,
                LDA,
                BI,
                LDA,
                MM,
                ALPHA,
                BETA,
                Q,
                LDA,
                Z,
                LDA,
                PL,
                DIFEST,
                WORK,
                LWORK,
                RWORK,
                IWORK,
                LIWORK,
                BWORK,
                LINFO);

            if (LINFO.value != 0 && LINFO.value != mn.MPLUSN + 2) {
              RESULT[1] = ULPINV;
              NOUT.println(
                  ' ZDRGSX: ZGGESX returned INFO=${LINFO.value.i6}.\n${' ' * 9}mn.N=${mn.MPLUSN.i6}, JTYPE=${PRTYPE.i6})');
              INFO.value = LINFO.value;
              continue;
            }

            // Compute the norm(A, B)

            zlacpy('Full', mn.MPLUSN, mn.MPLUSN, AI, LDA, WORK.asMatrix(),
                mn.MPLUSN);
            zlacpy('Full', mn.MPLUSN, mn.MPLUSN, BI, LDA,
                WORK(mn.MPLUSN * mn.MPLUSN + 1).asMatrix(), mn.MPLUSN);
            ABNRM = zlange('Fro', mn.MPLUSN, 2 * mn.MPLUSN, WORK.asMatrix(),
                mn.MPLUSN, RWORK);

            // Do tests (1) to (4)

            RESULT[2] = ZERO;
            zget51(1, mn.MPLUSN, A, LDA, AI, LDA, Q, LDA, Z, LDA, WORK, RWORK,
                RESULT(1));
            zget51(1, mn.MPLUSN, B, LDA, BI, LDA, Q, LDA, Z, LDA, WORK, RWORK,
                RESULT(2));
            zget51(3, mn.MPLUSN, B, LDA, BI, LDA, Q, LDA, Q, LDA, WORK, RWORK,
                RESULT(3));
            zget51(3, mn.MPLUSN, B, LDA, BI, LDA, Z, LDA, Z, LDA, WORK, RWORK,
                RESULT(4));
            NTEST = 4;

            // Do tests (5) and (6): check Schur form of A and
            // compare eigenvalues with diagonals.

            TEMP1 = ZERO;
            RESULT[5] = ZERO;
            RESULT[6] = ZERO;

            for (J = 1; J <= mn.MPLUSN; J++) {
              // 10
              ILABAD = false;
              TEMP2 = (ABS1(ALPHA[J] - AI[J][J]) /
                          max(SMLNUM, max(ABS1(ALPHA[J]), ABS1(AI[J][J]))) +
                      ABS1(BETA[J] - BI[J][J]) /
                          max(SMLNUM, max(ABS1(BETA[J]), ABS1(BI[J][J])))) /
                  ULP;
              if (J < mn.MPLUSN) {
                if (AI[J + 1][J] != Complex.zero) {
                  ILABAD = true;
                  RESULT[5] = ULPINV;
                }
              }
              if (J > 1) {
                if (AI[J][J - 1] != Complex.zero) {
                  ILABAD = true;
                  RESULT[5] = ULPINV;
                }
              }
              TEMP1 = max(TEMP1, TEMP2);
              if (ILABAD) {
                _print9997(NOUT, J, mn.MPLUSN, PRTYPE);
              }
            } // 10
            RESULT[6] = TEMP1;
            NTEST = NTEST + 2;

            // Test (7) (if sorting worked)

            RESULT[7] = ZERO;
            if (LINFO.value == mn.MPLUSN + 3) {
              RESULT[7] = ULPINV;
            } else if (MM.value != mn.N) {
              RESULT[7] = ULPINV;
            }
            NTEST++;

            // Test (8): compare the estimated value DIF and its
            // value. first, compute the exact DIF.

            RESULT[8] = ZERO;
            MN2 = MM.value * (mn.MPLUSN - MM.value) * 2;
            if (IFUNC >= 2 && MN2 <= NCMAX * NCMAX) {
              // Note: for either following two cases, there are
              // almost same number of test cases fail the test.

              zlakf2(
                  MM.value,
                  mn.MPLUSN - MM.value,
                  AI,
                  LDA,
                  AI(MM.value + 1, MM.value + 1),
                  BI,
                  BI(MM.value + 1, MM.value + 1),
                  C,
                  LDC);

              zgesvd('mn.N', 'mn.N', MN2, MN2, C, LDC, S, WORK.asMatrix(), 1,
                  WORK(2).asMatrix(), 1, WORK(3), LWORK - 2, RWORK, INFO);
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
            if (LINFO.value == (mn.MPLUSN + 2)) {
              if (DIFTRU > ABNRM * ULP) RESULT[9] = ULPINV;
              if ((IFUNC > 1) && (DIFEST[2] != ZERO)) RESULT[9] = ULPINV;
              if ((IFUNC == 1) && (PL[1] != ZERO)) RESULT[9] = ULPINV;
              NTEST++;
            }

            NTESTT = NTESTT + NTEST;

            // Print out tests which fail.

            for (J = 1; J <= 9; J++) {
              // 20
              if (RESULT[J] >= THRESH) {
                // If this is the first test to fail,
                // print a header to the data file.

                if (NERRS == 0) {
                  _print9996(NOUT, 'ZGX');

                  // Matrix types

                  NOUT.println(
                      ' Matrix types: \n  1:  A is a block diagonal matrix of Jordan blocks and B is the identity \n      matrix, \n  2:  A and B are upper triangular matrices, \n  3:  A and B are as type 2, but each second diagonal block in A_11 and \n      each third diagonal block in A_22 are 2x2 blocks,\n  4:  A and B are block diagonal matrices, \n  5:  (A,B) has potentially close or common eigenvalues.\n');

                  // Tests performed

                  _print9993(NOUT, 'unitary', 'transpose');
                }
                NERRS++;
                if (RESULT[J] < 10000.0) {
                  NOUT.println(
                      ' Matrix order=${mn.MPLUSN.i2}, type=${PRTYPE.i2}, a=${WEIGHT.d10_3}, order(A_11)=${mn.M.i2}, result ${J.i2} is ${RESULT[J].f8_2}');
                } else {
                  NOUT.println(
                      ' Matrix order=${mn.MPLUSN.i2}, type=${PRTYPE.i2}, a=${WEIGHT.d10_3}, order(A_11)=${mn.M.i2}, result ${J.i2} is ${RESULT[J].d10_3}');
                }
              }
            } // 20
          } // 30
        } // 40
      } // 50
    } // 60
  } else {
    // Read in data from file to check accuracy of condition estimation
    // Read input data until mn.N=0

    NPTKNT = 0;

    // } // 80
    try {
      while (true) {
        mn.MPLUSN = await NIN.readInt();
        if (mn.MPLUSN == 0) break;

        mn.N = await NIN.readInt();
        await NIN.readMatrix(AI, mn.MPLUSN, mn.MPLUSN);
        await NIN.readMatrix(BI, mn.MPLUSN, mn.MPLUSN);
        (PLTRU, DIFTRU) = await NIN.readDouble2();

        NPTKNT++;
        mn.FS = true;
        mn.K = 0;
        mn.M = mn.MPLUSN - mn.N;

        zlacpy('Full', mn.MPLUSN, mn.MPLUSN, AI, LDA, A, LDA);
        zlacpy('Full', mn.MPLUSN, mn.MPLUSN, BI, LDA, B, LDA);

        // Compute the Schur factorization while swapping the
        // m-by-m (1,1)-blocks with n-by-n (2,2)-blocks.

        zggesx(
            'V',
            'V',
            'S',
            zlctsx,
            'B',
            mn.MPLUSN,
            AI,
            LDA,
            BI,
            LDA,
            MM,
            ALPHA,
            BETA,
            Q,
            LDA,
            Z,
            LDA,
            PL,
            DIFEST,
            WORK,
            LWORK,
            RWORK,
            IWORK,
            LIWORK,
            BWORK,
            LINFO);

        if (LINFO.value != 0 && LINFO.value != mn.MPLUSN + 2) {
          RESULT[1] = ULPINV;
          NOUT.println(
              ' ZDRGSX: ZGGESX returned INFO=${LINFO.value.i6}.\n${' ' * 9}mn.N=${mn.MPLUSN.i6}, Input Example #${NPTKNT.i2})');
          continue;
        }

        // Compute the norm(A, B)
        //    (should this be norm of (A,B) or (AI,BI)?)

        zlacpy(
            'Full', mn.MPLUSN, mn.MPLUSN, AI, LDA, WORK.asMatrix(), mn.MPLUSN);
        zlacpy('Full', mn.MPLUSN, mn.MPLUSN, BI, LDA,
            WORK(mn.MPLUSN * mn.MPLUSN + 1).asMatrix(), mn.MPLUSN);
        ABNRM = zlange(
            'Fro', mn.MPLUSN, 2 * mn.MPLUSN, WORK.asMatrix(), mn.MPLUSN, RWORK);

        // Do tests (1) to (4)

        zget51(1, mn.MPLUSN, A, LDA, AI, LDA, Q, LDA, Z, LDA, WORK, RWORK,
            RESULT(1));
        zget51(1, mn.MPLUSN, B, LDA, BI, LDA, Q, LDA, Z, LDA, WORK, RWORK,
            RESULT(2));
        zget51(3, mn.MPLUSN, B, LDA, BI, LDA, Q, LDA, Q, LDA, WORK, RWORK,
            RESULT(3));
        zget51(3, mn.MPLUSN, B, LDA, BI, LDA, Z, LDA, Z, LDA, WORK, RWORK,
            RESULT(4));

        // Do tests (5) and (6): check Schur form of A and compare
        // eigenvalues with diagonals.

        NTEST = 6;
        TEMP1 = ZERO;
        RESULT[5] = ZERO;
        RESULT[6] = ZERO;

        for (J = 1; J <= mn.MPLUSN; J++) {
          // 110
          ILABAD = false;
          TEMP2 = (ABS1(ALPHA[J] - AI[J][J]) /
                      max(SMLNUM, max(ABS1(ALPHA[J]), ABS1(AI[J][J]))) +
                  ABS1(BETA[J] - BI[J][J]) /
                      max(SMLNUM, max(ABS1(BETA[J]), ABS1(BI[J][J])))) /
              ULP;
          if (J < mn.MPLUSN) {
            if (AI[J + 1][J] != Complex.zero) {
              ILABAD = true;
              RESULT[5] = ULPINV;
            }
          }
          if (J > 1) {
            if (AI[J][J - 1] != Complex.zero) {
              ILABAD = true;
              RESULT[5] = ULPINV;
            }
          }
          TEMP1 = max(TEMP1, TEMP2);
          if (ILABAD) {
            _print9997(NOUT, J, mn.MPLUSN, NPTKNT);
          }
        } // 110
        RESULT[6] = TEMP1;

        // Test (7) (if sorting worked)  <--------- need to be checked.

        NTEST = 7;
        RESULT[7] = ZERO;
        if (LINFO.value == mn.MPLUSN + 3) RESULT[7] = ULPINV;

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
        if (LINFO.value == (mn.MPLUSN + 2)) {
          if (DIFTRU > ABNRM * ULP) RESULT[9] = ULPINV;
          if ((IFUNC > 1) && (DIFEST[2] != ZERO)) RESULT[9] = ULPINV;
          if ((IFUNC == 1) && (PL[1] != ZERO)) RESULT[9] = ULPINV;
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

        NTESTT = NTESTT + NTEST;

        // Print out tests which fail.

        for (J = 1; J <= NTEST; J++) {
          // 120
          if (RESULT[J] >= THRESH) {
            // If this is the first test to fail,
            // print a header to the data file.

            if (NERRS == 0) {
              _print9996(NOUT, 'ZGX');

              // Matrix types

              NOUT.println('Input Example');

              // Tests performed

              _print9993(NOUT, 'unitary', 'transpose');
            }
            NERRS++;
            if (RESULT[J] < 10000.0) {
              NOUT.println(
                  ' Input example #${NPTKNT.i2}, matrix order=${mn.MPLUSN.i4}, result ${J.i2} is${RESULT[J].f8_2}');
            } else {
              NOUT.println(
                  ' Input example #${NPTKNT.i2}, matrix order=${mn.MPLUSN.i4}, result ${J.i2} is${(RESULT[J] * 10).d10_3}');
            }
          }
        } // 120

        // } // 130
      }
    } catch (_) {} // 140
  } // 150

  // Summary

  alasvm('ZGX', NOUT, NERRS, NTESTT, 0);

  WORK[1] = MAXWRK.toComplex();
}

void _print9996(Nout nout, String s) {
  nout.println(
      '\n ${s.a3} -- Complex Expert Generalized Schur form problem driver');
}

void _print9993(Nout nout, String s1, String s2) {
  nout.println(
      '\n Tests performed:  (S is Schur, T is triangular, Q and Z are $s1,\n${' ' * 19} a is alpha, b is beta, and \' means $s2.)\n  1 = | A - Q S Z\' | / ( |A| n ulp )      2 = | B - Q T Z\' | / ( |B| n ulp )\n  3 = | I - QQ\' | / ( n ulp )             4 = | I - ZZ\' | / ( n ulp )\n  5 = 1/ULP  if A is not in Schur form S\n  6 = difference between (alpha,beta) and diagonals of (S,T)\n  7 = 1/ULP  if SDIM is not the correct number of selected eigenvalues\n  8 = 1/ULP  if DIFEST/DIFTRU > 10*THRESH or DIFTRU/DIFEST > 10*THRESH\n  9 = 1/ULP  if DIFEST <> 0 or DIFTRU > ULP*norm(A,B) when reordering fails\n 10 = 1/ULP  if PLEST/PLTRU > THRESH or PLTRU/PLEST > THRESH\n    ( Test 10 is only for input examples )\n');
}

void _print9997(Nout nout, int i, int n, int jtype) {
  nout.println(
      ' ZDRGSX: S not in Schur form at eigenvalue ${i.i6}.\n${' ' * 9}mn.N=${n.i6}, JTYPE=${jtype.i6})');
}
