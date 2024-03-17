import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgges3.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlarfg.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dorm2r.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';

import '../matgen/dlarnd.dart';
import 'alasvm.dart';
import 'dget51.dart';
import 'dget53.dart';
import 'dget54.dart';
import 'dlatm4.dart';
import 'dlctes.dart';
import 'xlaenv.dart';

void ddrges3(
  final int NSIZES,
  final Array<int> NN_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final Array<int> ISEED_,
  final double THRESH,
  final Nout NOUNIT,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final Matrix<double> S_,
  final Matrix<double> T_,
  final Matrix<double> Q_,
  final int LDQ,
  final Matrix<double> Z_,
  final Array<double> ALPHAR_,
  final Array<double> ALPHAI_,
  final Array<double> BETA_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<double> RESULT_,
  final Array<bool> BWORK_,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final NN = NN_.having();
  final DOTYPE = DOTYPE_.having();
  final ISEED = ISEED_.having();
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDA);
  final S = S_.having(ld: LDA);
  final T = T_.having(ld: LDA);
  final Q = Q_.having(ld: LDQ);
  final Z = Z_.having(ld: LDQ);
  final ALPHAR = ALPHAR_.having();
  final ALPHAI = ALPHAI_.having();
  final BETA = BETA_.having();
  final WORK = WORK_.having();
  final RESULT = RESULT_.having();
  final BWORK = BWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  const MAXTYP = 26;
  bool BADNN, ILABAD;
  String SORT;
  int I,
      I1,
      IADD,
      IN,
      ISORT,
      J,
      JC,
      JR,
      JSIZE,
      JTYPE,
      KNTEIG,
      MAXWRK = 0,
      MINWRK,
      MTYPES,
      N,
      N1,
      NB,
      NERRS,
      NMAX,
      NTEST,
      NTESTT,
      RSUB;
  double SAFMAX, SAFMIN, TEMP1, ULP, ULPINV;
  final IERR = Box(0), IINFO = Box(0);
  final TEMP2 = Box(0.0);
  final IOLDSD = Array<int>(4), SDIM = Box(0);
  final RMAGN = Array<double>(4, offset: zeroIndexedArrayOffset);
  final KCLASS = Array.fromList([
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    3, //
  ]);
  final KZ1 = Array.fromList([0, 1, 2, 1, 3, 3]);
  final KZ2 = Array.fromList([0, 0, 1, 2, 1, 1]);
  final KADD = Array.fromList([0, 0, 0, 0, 3, 2]);
  final KATYPE = Array.fromList([
    0, 1, 0, 1, 2, 3, 4, 1, 4, 4, 1, 1, 4, 4, 4, 2, 4, 5, 8, 7, 9, 4, 4, 4, 4,
    0, //
  ]);
  final KBTYPE = Array.fromList([
    0, 0, 1, 1, 2, -3, 1, 4, 1, 1, 4, 4, 1, 1, -4, 2, -4, 8, 8, 8, 8, 8, 8, 8,
    8, 0, //
  ]);
  final KAZERO = Array.fromList([
    1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 1, 2, 2, 3, 1, 3, 5, 5, 5, 5, 3, 3, 3, 3,
    1, //
  ]);
  final KBZERO = Array.fromList([
    1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 4, 1, 4, 6, 6, 6, 6, 4, 4, 4, 4,
    1, //
  ]);
  final KAMAGN = Array.fromList([
    1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 2, 3, 2, 3, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 2,
    1, //
  ]);
  final KBMAGN = Array.fromList([
    1, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2, 2, 3, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2,
    1, //
  ]);
  final KTRIAN = Array.fromList([
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, //
  ]);
  final IASIGN = Array.fromList([
    0, 0, 0, 0, 0, 0, 2, 0, 2, 2, 0, 0, 2, 2, 2, 0, 2, 0, 0, 0, 2, 2, 2, 2, 2,
    0, //
  ]);
  final IBSIGN = Array.fromList([
    0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 2, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    0, //
  ]);

  // Check for errors

  INFO.value = 0;

  BADNN = false;
  NMAX = 1;
  for (J = 1; J <= NSIZES; J++) {
    NMAX = max(NMAX, NN[J]);
    if (NN[J] < 0) BADNN = true;
  }

  if (NSIZES < 0) {
    INFO.value = -1;
  } else if (BADNN) {
    INFO.value = -2;
  } else if (NTYPES < 0) {
    INFO.value = -3;
  } else if (THRESH < ZERO) {
    INFO.value = -6;
  } else if (LDA <= 1 || LDA < NMAX) {
    INFO.value = -9;
  } else if (LDQ <= 1 || LDQ < NMAX) {
    INFO.value = -14;
  }

  // Compute workspace
  // (Note: Comments in the code beginning "Workspace:" describe the
  // minimal amount of workspace needed at that point in the code,
  // as well as the preferred amount for good performance.
  // NB refers to the optimal block size for the immediately
  // following subroutine, as returned by ILAENV.

  MINWRK = 1;
  if (INFO.value == 0 && LWORK >= 1) {
    MINWRK = max(10 * (NMAX + 1), 3 * NMAX * NMAX);
    NB = max(
        max(1, ilaenv(1, 'DGEQRF', ' ', NMAX, NMAX, -1, -1)),
        max(
          ilaenv(1, 'DORMQR', 'LT', NMAX, NMAX, NMAX, -1),
          ilaenv(1, 'DORGQR', ' ', NMAX, NMAX, NMAX, -1),
        ));
    MAXWRK = max(10 * (NMAX + 1), max(2 * NMAX + NMAX * NB, 3 * NMAX * NMAX));
    WORK[1] = MAXWRK.toDouble();
  }

  if (LWORK < MINWRK) INFO.value = -20;

  if (INFO.value != 0) {
    xerbla('DDRGES3', -INFO.value);
    return;
  }

  // Quick return if possible

  if (NSIZES == 0 || NTYPES == 0) return;

  SAFMIN = dlamch('Safe minimum');
  ULP = dlamch('Epsilon') * dlamch('Base');
  SAFMIN /= ULP;
  SAFMAX = ONE / SAFMIN;
  ULPINV = ONE / ULP;

  // The values RMAGN[2:3] depend on N, see below.

  RMAGN[0] = ZERO;
  RMAGN[1] = ONE;

  // Loop over matrix sizes

  NTESTT = 0;
  NERRS = 0;

  for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) {
    N = NN[JSIZE];
    N1 = max(1, N);
    RMAGN[2] = SAFMAX * ULP / N1.toDouble();
    RMAGN[3] = SAFMIN * ULPINV * N1.toDouble();

    if (NSIZES != 1) {
      MTYPES = min(MAXTYP, NTYPES);
    } else {
      MTYPES = min(MAXTYP + 1, NTYPES);
    }

    // Loop over matrix types

    for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
      if (!DOTYPE[JTYPE]) continue;
      NTEST = 0;

      // Save ISEED in case of an error.

      for (J = 1; J <= 4; J++) {
        IOLDSD[J] = ISEED[J];
      }

      // Initialize RESULT

      for (J = 1; J <= 13; J++) {
        RESULT[J] = ZERO;
      }

      // Generate test matrices A and B

      // Description of control parameters:
      //
      // KZLASS: =1 means w/o rotation, =2 means w/ rotation,
      //         =3 means random.
      // KATYPE: the "type" to be passed to DLATM4 for computing A.
      // KAZERO: the pattern of zeros on the diagonal for A:
      //         =1: ( xxx ), =2: (0, xxx ) =3: ( 0, 0, xxx, 0 ),
      //         =4: ( 0, xxx, 0, 0 ), =5: ( 0, 0, 1, xxx, 0 ),
      //         =6: ( 0, 1, 0, xxx, 0 ).  (xxx means a string of
      //         non-zero entries.)
      // KAMAGN: the magnitude of the matrix: =0: zero, =1: O(1),
      //         =2: large, =3: small.
      // IASIGN: 1 if the diagonal elements of A are to be
      //         multiplied by a random magnitude 1 number, =2 if
      //         randomly chosen diagonal blocks are to be rotated
      //         to form 2x2 blocks.
      // KBTYPE, KBZERO, KBMAGN, IBSIGN: the same, but for B.
      // KTRIAN: =0: don't fill in the upper triangle, =1: do.
      // KZ1, KZ2, KADD: used to implement KAZERO and KBZERO.
      // RMAGN: used to implement KAMAGN and KBMAGN.

      if (MTYPES <= MAXTYP) {
        IINFO.value = 0;
        if (KCLASS[JTYPE] < 3) {
          // Generate A (w/o rotation)

          if (KATYPE[JTYPE].abs() == 3) {
            IN = 2 * ((N - 1) ~/ 2) + 1;
            if (IN != N) dlaset('Full', N, N, ZERO, ZERO, A, LDA);
          } else {
            IN = N;
          }
          dlatm4(
              KATYPE[JTYPE],
              IN,
              KZ1[KAZERO[JTYPE]],
              KZ2[KAZERO[JTYPE]],
              IASIGN[JTYPE],
              RMAGN[KAMAGN[JTYPE]],
              ULP,
              RMAGN[KTRIAN[JTYPE] * KAMAGN[JTYPE]],
              2,
              ISEED,
              A,
              LDA);
          IADD = KADD[KAZERO[JTYPE]];
          if (IADD > 0 && IADD <= N) A[IADD][IADD] = ONE;

          // Generate B (w/o rotation)

          if (KBTYPE[JTYPE].abs() == 3) {
            IN = 2 * ((N - 1) ~/ 2) + 1;
            if (IN != N) dlaset('Full', N, N, ZERO, ZERO, B, LDA);
          } else {
            IN = N;
          }
          dlatm4(
              KBTYPE[JTYPE],
              IN,
              KZ1[KBZERO[JTYPE]],
              KZ2[KBZERO[JTYPE]],
              IBSIGN[JTYPE],
              RMAGN[KBMAGN[JTYPE]],
              ONE,
              RMAGN[KTRIAN[JTYPE] * KBMAGN[JTYPE]],
              2,
              ISEED,
              B,
              LDA);
          IADD = KADD[KBZERO[JTYPE]];
          if (IADD != 0 && IADD <= N) B[IADD][IADD] = ONE;

          if (KCLASS[JTYPE] == 2 && N > 0) {
            // Include rotations

            // Generate Q, Z as Householder transformations times
            // a diagonal matrix.

            for (JC = 1; JC <= N - 1; JC++) {
              for (JR = JC; JR <= N; JR++) {
                Q[JR][JC] = dlarnd(3, ISEED);
                Z[JR][JC] = dlarnd(3, ISEED);
              }
              dlarfg(N + 1 - JC, Q.box(JC, JC), Q(JC + 1, JC).asArray(), 1,
                  WORK.box(JC));
              WORK[2 * N + JC] = sign(ONE, Q[JC][JC]).toDouble();
              Q[JC][JC] = ONE;
              dlarfg(N + 1 - JC, Z.box(JC, JC), Z(JC + 1, JC).asArray(), 1,
                  WORK.box(N + JC));
              WORK[3 * N + JC] = sign(ONE, Z[JC][JC]).toDouble();
              Z[JC][JC] = ONE;
            }
            Q[N][N] = ONE;
            WORK[N] = ZERO;
            WORK[3 * N] = sign(ONE, dlarnd(2, ISEED)).toDouble();
            Z[N][N] = ONE;
            WORK[2 * N] = ZERO;
            WORK[4 * N] = sign(ONE, dlarnd(2, ISEED)).toDouble();

            // Apply the diagonal matrices

            for (JC = 1; JC <= N; JC++) {
              for (JR = 1; JR <= N; JR++) {
                A[JR][JC] = WORK[2 * N + JR] * WORK[3 * N + JC] * A[JR][JC];
                B[JR][JC] = WORK[2 * N + JR] * WORK[3 * N + JC] * B[JR][JC];
              }
            }
            dorm2r('L', 'N', N, N, N - 1, Q, LDQ, WORK, A, LDA, WORK(2 * N + 1),
                IINFO);
            if (IINFO.value == 0) {
              dorm2r('R', 'T', N, N, N - 1, Z, LDQ, WORK(N + 1), A, LDA,
                  WORK(2 * N + 1), IINFO);
              if (IINFO.value == 0) {
                dorm2r('L', 'N', N, N, N - 1, Q, LDQ, WORK, B, LDA,
                    WORK(2 * N + 1), IINFO);
                if (IINFO.value == 0) {
                  dorm2r('R', 'T', N, N, N - 1, Z, LDQ, WORK(N + 1), B, LDA,
                      WORK(2 * N + 1), IINFO);
                }
              }
            }
          }
        } else {
          // Random matrices

          for (JC = 1; JC <= N; JC++) {
            for (JR = 1; JR <= N; JR++) {
              A[JR][JC] = RMAGN[KAMAGN[JTYPE]] * dlarnd(2, ISEED);
              B[JR][JC] = RMAGN[KBMAGN[JTYPE]] * dlarnd(2, ISEED);
            }
          }
        }

        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'Generator', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }
      }

      for (I = 1; I <= 13; I++) {
        RESULT[I] = -ONE;
      }

      // Test with and without sorting of eigenvalues

      for (ISORT = 0; ISORT <= 1; ISORT++) {
        if (ISORT == 0) {
          SORT = 'N';
          RSUB = 0;
        } else {
          SORT = 'S';
          RSUB = 5;
        }

        // Call XLAENV to set the parameters used in DLAQZ0

        xlaenv(12, 10);
        xlaenv(13, 12);
        xlaenv(14, 13);
        xlaenv(15, 2);
        xlaenv(17, 10);

        // Call DGGES3 to compute S, T, Q, Z, alpha, and beta.

        dlacpy('Full', N, N, A, LDA, S, LDA);
        dlacpy('Full', N, N, B, LDA, T, LDA);
        NTEST = 1 + RSUB + ISORT;
        RESULT[1 + RSUB + ISORT] = ULPINV;
        dgges3('V', 'V', SORT, dlctes, N, S, LDA, T, LDA, SDIM, ALPHAR, ALPHAI,
            BETA, Q, LDQ, Z, LDQ, WORK, LWORK, BWORK, IINFO);
        if (IINFO.value != 0 && IINFO.value != N + 2) {
          RESULT[1 + RSUB + ISORT] = ULPINV;
          _print9999(NOUNIT, 'DGGES3', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break;
        }

        NTEST = 4 + RSUB;

        // Do tests 1--4 (or tests 7--9 when reordering )

        if (ISORT == 0) {
          dget51(1, N, A, LDA, S, LDA, Q, LDQ, Z, LDQ, WORK, RESULT.box(1));
          dget51(1, N, B, LDA, T, LDA, Q, LDQ, Z, LDQ, WORK, RESULT.box(2));
        } else {
          dget54(N, A, LDA, B, LDA, S, LDA, T, LDA, Q, LDQ, Z, LDQ, WORK,
              RESULT.box(7));
        }
        dget51(
            3, N, A, LDA, T, LDA, Q, LDQ, Q, LDQ, WORK, RESULT.box(3 + RSUB));
        dget51(
            3, N, B, LDA, T, LDA, Z, LDQ, Z, LDQ, WORK, RESULT.box(4 + RSUB));

        // Do test 5 and 6 (or Tests 10 and 11 when reordering):
        // check Schur form of A and compare eigenvalues with
        // diagonals.

        NTEST = 6 + RSUB;
        TEMP1 = ZERO;

        for (J = 1; J <= N; J++) {
          ILABAD = false;
          if (ALPHAI[J] == ZERO) {
            TEMP2.value = ((ALPHAR[J] - S[J][J]).abs() /
                        max(
                          SAFMIN,
                          max(ALPHAR[J].abs(), S[J][J].abs()),
                        ) +
                    (BETA[J] - T[J][J]).abs() /
                        max(SAFMIN, max(BETA[J].abs(), T[J][J].abs()))) /
                ULP;

            if (J < N) {
              if (S[J + 1][J] != ZERO) {
                ILABAD = true;
                RESULT[5 + RSUB] = ULPINV;
              }
            }
            if (J > 1) {
              if (S[J][J - 1] != ZERO) {
                ILABAD = true;
                RESULT[5 + RSUB] = ULPINV;
              }
            }
          } else {
            if (ALPHAI[J] > ZERO) {
              I1 = J;
            } else {
              I1 = J - 1;
            }
            if (I1 <= 0 || I1 >= N) {
              ILABAD = true;
            } else if (I1 < N - 1) {
              if (S[I1 + 2][I1 + 1] != ZERO) {
                ILABAD = true;
                RESULT[5 + RSUB] = ULPINV;
              }
            } else if (I1 > 1) {
              if (S[I1][I1 - 1] != ZERO) {
                ILABAD = true;
                RESULT[5 + RSUB] = ULPINV;
              }
            }
            if (!ILABAD) {
              dget53(S(I1, I1), LDA, T(I1, I1), LDA, BETA[J], ALPHAR[J],
                  ALPHAI[J], TEMP2, IERR);
              if (IERR.value >= 3) {
                NOUNIT.println(
                    ' DDRGES3: DGET53 returned INFOe=${IERR.value.i1} for eigenvalue ${J.i6}.\n${' ' * 9}N=${N.i6}, JTYPE=${JTYPE.i6}, ISEED=(${IOLDSD.i4(4, ',')})');
                INFO.value = (IERR.value).abs();
              }
            } else {
              TEMP2.value = ULPINV;
            }
          }
          TEMP1 = max(TEMP1, TEMP2.value);
          if (ILABAD) {
            NOUNIT.println(
                ' DDRGES3: S not in Schur form at eigenvalue ${J.i6}.\n${' ' * 9}N=${N.i6}, JTYPE=${JTYPE.i6}, ISEED=(${IOLDSD.i5(4, ',')})');
          }
        }
        RESULT[6 + RSUB] = TEMP1;

        if (ISORT >= 1) {
          // Do test 12

          NTEST = 12;
          RESULT[12] = ZERO;
          KNTEIG = 0;
          for (I = 1; I <= N; I++) {
            if (dlctes(ALPHAR[I], ALPHAI[I], BETA[I]) ||
                dlctes(ALPHAR[I], -ALPHAI[I], BETA[I])) {
              KNTEIG++;
            }
            if (I < N) {
              if ((dlctes(ALPHAR[I + 1], ALPHAI[I + 1], BETA[I + 1]) ||
                      dlctes(ALPHAR[I + 1], -ALPHAI[I + 1], BETA[I + 1])) &&
                  (!(dlctes(ALPHAR[I], ALPHAI[I], BETA[I]) ||
                      dlctes(ALPHAR[I], -ALPHAI[I], BETA[I]))) &&
                  IINFO.value != N + 2) {
                RESULT[12] = ULPINV;
              }
            }
          }
          if (SDIM.value != KNTEIG) {
            RESULT[12] = ULPINV;
          }
        }
      }

      // End of Loop -- Check for RESULT[j] > THRESH

      NTESTT += NTEST;

      // Print out tests which fail.

      for (JR = 1; JR <= NTEST; JR++) {
        if (RESULT[JR] >= THRESH) {
          // If this is the first test to fail,
          // print a header to the data file.

          if (NERRS == 0) {
            NOUNIT.println('\n DGS -- Real Generalized Schur form driver');

            // Matrix types

            NOUNIT.println(' Matrix types (see DDRGES3 for details): ');
            NOUNIT.println(
                ' Special Matrices:${' ' * 23}(J\'=transposed Jordan block)\n   1=(0,0)  2=(I,0)  3=(0,I)  4=(I,I)  5=(J\',J\')  6=(diag(J\',I), diag(I,J\'))\n Diagonal Matrices:  ( D=diag(0,1,2,...) )\n   7=(D,I)   9=(large*D, small*I)  11=(large*I, small*D)  13=(large*D, large*I)\n   8=(I,D)  10=(small*D, large*I)  12=(small*I, large*D)  14=(small*D, small*I)\n  15=(D, reversed D)');
            NOUNIT.println(
                ' Matrices Rotated by Random Orthogonal Matrices U, V:\n  16=Transposed Jordan Blocks             19=geometric alpha, beta=0,1\n  17=arithm. alpha&beta                   20=arithmetic alpha, beta=0,1\n  18=clustered alpha, beta=0,1            21=random alpha, beta=0,1\n Large & Small Matrices:\n  22=(large, small)   23=(small,large)    24=(small,small)    25=(large,large)\n  26=random O(1) matrices.');

            // Tests performed

            NOUNIT.println(
                '\n Tests performed:  (S is Schur, T is triangular, Q and Z are orthogonal,\n${' ' * 19}l and r are the appropriate left and right\n${' ' * 19}eigenvectors, resp., a is alpha, b is beta, and\n${' ' * 19}\' means transpose.)\n Without ordering: \n  1 = | A - Q S Z\' | / ( |A| n ulp )      2 = | B - Q T Z\' | / ( |B| n ulp )\n  3 = | I - QQ\' | / ( n ulp )             4 = | I - ZZ\' | / ( n ulp )\n  5 = A is in Schur form S\n  6 = difference between (alpha,beta) and diagonals of (S,T)\n With ordering: \n  7 = | (A,B) - Q (S,T) Z\' | / ( |(A,B)| n ulp )  \n  8 = | I - QQ\' | / ( n ulp )            9 = | I - ZZ\' | / ( n ulp )\n 10 = A is in Schur form S\n 11 = difference between (alpha,beta) and diagonals of (S,T)\n 12 = SDIM is the correct number of selected eigenvalues\n');
          }
          NERRS++;
          if (RESULT[JR] < 10000.0) {
            NOUNIT.println(
                ' Matrix order=${N.i5}, type=${JTYPE.i2}, seed=${IOLDSD.i4(4, ',')} result ${JR.i2} is${RESULT[JR].f8_2}');
          } else {
            NOUNIT.println(
                ' Matrix order=${N.i5}, type=${JTYPE.i2}, seed=${IOLDSD.i4(4, ',')} result ${JR.i2} is${(RESULT[JR] * 10).d10_3}');
          }
        }
      }
    }
  }

  // Summary
  alasvm('DGS', NOUNIT, NERRS, NTESTT, 0);

  WORK[1] = MAXWRK.toDouble();
}

void _print9999(
  final Nout nout,
  final String s,
  final int info,
  final int n,
  final int jtype,
  final Array<int> iseed,
) {
  nout.println(
      ' DDRGES3: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i4(4, ',')})');
}
