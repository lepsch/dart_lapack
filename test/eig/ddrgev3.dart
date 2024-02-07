import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dggev3.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlarfg.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dorm2r.dart';
import 'package:lapack/src/f2c/sign.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';

import '../matgen/dlarnd.dart';
import 'alasvm.dart';
import 'dget52.dart';
import 'dlatm4.dart';
import 'xlaenv.dart';

void ddrgev3(
  final int NSIZES,
  final Array<int> NN,
  final int NTYPES,
  final Array<bool> DOTYPE,
  final Array<int> ISEED,
  final double THRESH,
  final Nout NOUNIT,
  final Matrix<double> A,
  final int LDA,
  final Matrix<double> B,
  final Matrix<double> S,
  final Matrix<double> T,
  final Matrix<double> Q,
  final int LDQ,
  final Matrix<double> Z,
  final Matrix<double> QE,
  final int LDQE,
  final Array<double> ALPHAR,
  final Array<double> ALPHAI,
  final Array<double> BETA,
  final Array<double> ALPHR1,
  final Array<double> ALPHI1,
  final Array<double> BETA1,
  final Array<double> WORK,
  final int LWORK,
  final Array<double> RESULT,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  const MAXTYP = 27;
  bool BADNN;
  int I,
      IADD,
      IN,
      J,
      JC,
      JR,
      JSIZE,
      JTYPE,
      MAXWRK = 0,
      MINWRK,
      MTYPES,
      N,
      N1,
      NERRS,
      NMATS,
      NMAX,
      NTESTT;
  double SAFMAX, SAFMIN, ULP, ULPINV;
  final IERR = Box(0);
  final IOLDSD = Array<int>(4);
  final RMAGN = Array<double>(4, offset: 1);
  const KCLASS = [
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    3, 4, //
  ];
  const KZ1 = [0, 1, 2, 1, 3, 3];
  const KZ2 = [0, 0, 1, 2, 1, 1];
  const KADD = [0, 0, 0, 0, 3, 2];
  const KATYPE = [
    0, 1, 0, 1, 2, 3, 4, 1, 4, 4, 1, 1, 4, 4, 4, 2, 4, 5, 8, 7, 9, 4, 4, 4, 4,
    0, 0, //
  ];
  const KBTYPE = [
    0, 0, 1, 1, 2, -3, 1, 4, 1, 1, 4, 4, 1, 1, -4, 2, -4, 8, 8, 8, 8, 8, 8, 8,
    8, 0, 0, //
  ];
  const KAZERO = [
    1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 1, 2, 2, 3, 1, 3, 5, 5, 5, 5, 3, 3, 3, 3,
    1, 1, //
  ];
  const KBZERO = [
    1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 4, 1, 4, 6, 6, 6, 6, 4, 4, 4, 4,
    1, 1, //
  ];
  const KAMAGN = [
    1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 2, 3, 2, 3, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 2,
    1, 3, //
  ];
  const KBMAGN = [
    1, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2, 2, 3, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2,
    1, 3, //
  ];
  const KTRIAN = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, //
  ];
  const IASIGN = [
    0, 0, 0, 0, 0, 0, 2, 0, 2, 2, 0, 0, 2, 2, 2, 0, 2, 0, 0, 0, 2, 2, 2, 2, 2,
    0, 0, //
  ];
  const IBSIGN = [
    0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 2, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, //
  ];

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
  } else if (LDQE <= 1 || LDQE < NMAX) {
    INFO.value = -17;
  }

  // Compute workspace
  // (Note: Comments in the code beginning "Workspace:" describe the
  // minimal amount of workspace needed at that point in the code,
  // as well as the preferred amount for good performance.
  // NB refers to the optimal block size for the immediately
  // following subroutine, as returned by ILAENV.

  MINWRK = 1;
  if (INFO.value == 0 && LWORK >= 1) {
    MINWRK = max(1, max(8 * NMAX, NMAX * (NMAX + 1)));
    MAXWRK = 7 * NMAX + NMAX * ilaenv(1, 'DGEQRF', ' ', NMAX, 1, NMAX, 0);
    MAXWRK = max(MAXWRK, NMAX * (NMAX + 1));
    WORK[1] = MAXWRK.toDouble();
  }

  if (LWORK < MINWRK) INFO.value = -25;

  if (INFO.value != 0) {
    xerbla('DDRGEV3', -INFO.value);
    return;
  }

  // Quick return if possible

  if (NSIZES == 0 || NTYPES == 0) return;

  SAFMIN = dlamch('Safe minimum');
  ULP = dlamch('Epsilon') * dlamch('Base');
  SAFMIN = SAFMIN / ULP;
  SAFMAX = ONE / SAFMIN;
  ULPINV = ONE / ULP;

  // The values RMAGN(2:3) depend on N, see below.

  RMAGN[0] = ZERO;
  RMAGN[1] = ONE;

  // Loop over sizes, types

  NTESTT = 0;
  NERRS = 0;
  NMATS = 0;

  for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) {
    N = NN[JSIZE];
    N1 = max(1, N);
    RMAGN[2] = SAFMAX * ULP / N1.toDouble();
    RMAGN[3] = SAFMIN * ULPINV * N1;

    if (NSIZES != 1) {
      MTYPES = min(MAXTYP, NTYPES);
    } else {
      MTYPES = min(MAXTYP + 1, NTYPES);
    }

    for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
      if (!DOTYPE[JTYPE]) continue;
      NMATS = NMATS + 1;

      // Save ISEED in case of an error.

      for (J = 1; J <= 4; J++) {
        IOLDSD[J] = ISEED[J];
      }

      // Generate test matrices A and B

      // Description of control parameters:
      //
      // KZLASS: =1 means w/o rotation, =2 means w/ rotation,
      // =3 means random, =4 means random generalized
      // upper Hessenberg.
      // KATYPE: the "type" to be passed to DLATM4 for computing A.
      // KAZERO: the pattern of zeros on the diagonal for A:
      // =1: ( xxx ), =2: (0, xxx ) =3: ( 0, 0, xxx, 0 ),
      // =4: ( 0, xxx, 0, 0 ), =5: ( 0, 0, 1, xxx, 0 ),
      // =6: ( 0, 1, 0, xxx, 0 ).  (xxx means a string of
      // non-zero entries.)
      // KAMAGN: the magnitude of the matrix: =0: zero, =1: O(1),
      // =2: large, =3: small.
      // IASIGN: 1 if the diagonal elements of A are to be
      // multiplied by a random magnitude 1 number, =2 if
      // randomly chosen diagonal blocks are to be rotated
      // to form 2x2 blocks.
      // KBTYPE, KBZERO, KBMAGN, IBSIGN: the same, but for B.
      // KTRIAN: =0: don't fill in the upper triangle, =1: do.
      // KZ1, KZ2, KADD: used to implement KAZERO and KBZERO.
      // RMAGN: used to implement KAMAGN and KBMAGN.

      if (MTYPES <= MAXTYP) {
        IERR.value = 0;
        if (KCLASS[JTYPE] < 3) {
          // Generate A (w/o rotation)

          if ((KATYPE[JTYPE]).abs() == 3) {
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
            LDA,
          );
          IADD = KADD[KAZERO[JTYPE]];
          if (IADD > 0 && IADD <= N) A[IADD][IADD] = ONE;

          // Generate B (w/o rotation)

          if ((KBTYPE[JTYPE]).abs() == 3) {
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
            LDA,
          );
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
              dlarfg(
                N + 1 - JC,
                Q.box(JC, JC),
                Q(JC + 1, JC).asArray(),
                1,
                WORK.box(JC),
              );
              WORK[2 * N + JC] = sign(ONE, Q[JC][JC]).toDouble();
              Q[JC][JC] = ONE;
              dlarfg(
                N + 1 - JC,
                Z.box(JC, JC),
                Z(JC + 1, JC).asArray(),
                1,
                WORK.box(N + JC),
              );
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
            dorm2r(
              'L',
              'N',
              N,
              N,
              N - 1,
              Q,
              LDQ,
              WORK,
              A,
              LDA,
              WORK(2 * N + 1),
              IERR.value,
            );
            if (IERR.value == 0) {
              dorm2r(
                'R',
                'T',
                N,
                N,
                N - 1,
                Z,
                LDQ,
                WORK(N + 1),
                A,
                LDA,
                WORK(2 * N + 1),
                IERR.value,
              );
              if (IERR.value == 0) {
                dorm2r(
                  'L',
                  'N',
                  N,
                  N,
                  N - 1,
                  Q,
                  LDQ,
                  WORK,
                  B,
                  LDA,
                  WORK(2 * N + 1),
                  IERR.value,
                );
                if (IERR.value == 0) {
                  dorm2r(
                    'R',
                    'T',
                    N,
                    N,
                    N - 1,
                    Z,
                    LDQ,
                    WORK(N + 1),
                    B,
                    LDA,
                    WORK(2 * N + 1),
                    IERR.value,
                  );
                }
              }
            }
          }
        } else if (KCLASS[JTYPE] == 3) {
          // Random matrices

          for (JC = 1; JC <= N; JC++) {
            for (JR = 1; JR <= N; JR++) {
              A[JR][JC] = RMAGN[KAMAGN[JTYPE]] * dlarnd(2, ISEED);
              B[JR][JC] = RMAGN[KBMAGN[JTYPE]] * dlarnd(2, ISEED);
            }
          }
        } else {
          // Random upper Hessenberg pencil with singular B

          for (JC = 1; JC <= N; JC++) {
            for (JR = 1; JR <= min(JC + 1, N); JR++) {
              A[JR][JC] = RMAGN[KAMAGN[JTYPE]] * dlarnd(2, ISEED);
            }
            for (JR = JC + 2; JR <= N; JR++) {
              A[JR][JC] = ZERO;
            }
          }
          for (JC = 1; JC <= N; JC++) {
            for (JR = 1; JR <= JC; JR++) {
              B[JR][JC] = RMAGN[KAMAGN[JTYPE]] * dlarnd(2, ISEED);
            }
            for (JR = JC + 1; JR <= N; JR++) {
              B[JR][JC] = ZERO;
            }
          }
          for (JC = 1; JC <= N; JC += 4) {
            B[JC][JC] = ZERO;
          }
        }

        if (IERR.value != 0) {
          _print9999(NOUNIT, 'Generator', IERR.value, N, JTYPE, IOLDSD);
          INFO.value = (IERR.value).abs();
          return;
        }
      }

      for (I = 1; I <= 7; I++) {
        RESULT[I] = -ONE;
      }

      // Call XLAENV to set the parameters used in DLAQZ0

      xlaenv(12, 10);
      xlaenv(13, 12);
      xlaenv(14, 13);
      xlaenv(15, 2);
      xlaenv(17, 10);

      // Call DGGEV3 to compute eigenvalues and eigenvectors.
      while (true) {
        dlacpy(' ', N, N, A, LDA, S, LDA);
        dlacpy(' ', N, N, B, LDA, T, LDA);
        dggev3(
          'V',
          'V',
          N,
          S,
          LDA,
          T,
          LDA,
          ALPHAR,
          ALPHAI,
          BETA,
          Q,
          LDQ,
          Z,
          LDQ,
          WORK,
          LWORK,
          IERR,
        );
        if (IERR.value != 0 && IERR.value != N + 1) {
          RESULT[1] = ULPINV;
          _print9999(NOUNIT, 'DGGEV31', IERR.value, N, JTYPE, IOLDSD);
          INFO.value = (IERR.value).abs();
          break;
        }

        // Do the tests (1) and (2)

        dget52(
          true,
          N,
          A,
          LDA,
          B,
          LDA,
          Q,
          LDQ,
          ALPHAR,
          ALPHAI,
          BETA,
          WORK,
          RESULT(1),
        );
        if (RESULT[2] > THRESH) {
          _print9998(NOUNIT, 'Left', 'DGGEV31', RESULT[2], N, JTYPE, IOLDSD);
        }

        // Do the tests (3) and (4)

        dget52(
          false,
          N,
          A,
          LDA,
          B,
          LDA,
          Z,
          LDQ,
          ALPHAR,
          ALPHAI,
          BETA,
          WORK,
          RESULT(3),
        );
        if (RESULT[4] > THRESH) {
          _print9998(NOUNIT, 'Right', 'DGGEV31', RESULT[4], N, JTYPE, IOLDSD);
        }

        // Do the test (5)

        dlacpy(' ', N, N, A, LDA, S, LDA);
        dlacpy(' ', N, N, B, LDA, T, LDA);
        dggev3(
          'N',
          'N',
          N,
          S,
          LDA,
          T,
          LDA,
          ALPHR1,
          ALPHI1,
          BETA1,
          Q,
          LDQ,
          Z,
          LDQ,
          WORK,
          LWORK,
          IERR,
        );
        if (IERR.value != 0 && IERR.value != N + 1) {
          RESULT[1] = ULPINV;
          _print9999(NOUNIT, 'DGGEV32', IERR.value, N, JTYPE, IOLDSD);
          INFO.value = (IERR.value).abs();
          break;
        }

        for (J = 1; J <= N; J++) {
          if (ALPHAR[J] != ALPHR1[J] ||
              ALPHAI[J] != ALPHI1[J] ||
              BETA[J] != BETA1[J]) RESULT[5] = ULPINV;
        }

        // Do the test (6): Compute eigenvalues and left eigenvectors,
        // and test them

        dlacpy(' ', N, N, A, LDA, S, LDA);
        dlacpy(' ', N, N, B, LDA, T, LDA);
        dggev3(
          'V',
          'N',
          N,
          S,
          LDA,
          T,
          LDA,
          ALPHR1,
          ALPHI1,
          BETA1,
          QE,
          LDQE,
          Z,
          LDQ,
          WORK,
          LWORK,
          IERR,
        );
        if (IERR.value != 0 && IERR.value != N + 1) {
          RESULT[1] = ULPINV;
          _print9999(NOUNIT, 'DGGEV33', IERR.value, N, JTYPE, IOLDSD);
          INFO.value = (IERR.value).abs();
          break;
        }

        for (J = 1; J <= N; J++) {
          if (ALPHAR[J] != ALPHR1[J] ||
              ALPHAI[J] != ALPHI1[J] ||
              BETA[J] != BETA1[J]) RESULT[6] = ULPINV;
        }

        for (J = 1; J <= N; J++) {
          for (JC = 1; JC <= N; JC++) {
            if (Q[J][JC] != QE[J][JC]) RESULT[6] = ULPINV;
          }
        }

        // DO the test (7): Compute eigenvalues and right eigenvectors,
        // and test them

        dlacpy(' ', N, N, A, LDA, S, LDA);
        dlacpy(' ', N, N, B, LDA, T, LDA);
        dggev3(
          'N',
          'V',
          N,
          S,
          LDA,
          T,
          LDA,
          ALPHR1,
          ALPHI1,
          BETA1,
          Q,
          LDQ,
          QE,
          LDQE,
          WORK,
          LWORK,
          IERR,
        );
        if (IERR.value != 0 && IERR.value != N + 1) {
          RESULT[1] = ULPINV;
          _print9999(NOUNIT, 'DGGEV34', IERR.value, N, JTYPE, IOLDSD);
          INFO.value = (IERR.value).abs();
          break;
        }

        for (J = 1; J <= N; J++) {
          if (ALPHAR[J] != ALPHR1[J] ||
              ALPHAI[J] != ALPHI1[J] ||
              BETA[J] != BETA1[J]) RESULT[7] = ULPINV;
        }

        for (J = 1; J <= N; J++) {
          for (JC = 1; JC <= N; JC++) {
            if (Z[J][JC] != QE[J][JC]) RESULT[7] = ULPINV;
          }
        }

        // End of Loop -- Check for RESULT[j] > THRESH
        break;
      }

      NTESTT = NTESTT + 7;

      // Print out tests which fail.

      for (JR = 1; JR <= 7; JR++) {
        if (RESULT[JR] >= THRESH) {
          // If this is the first test to fail,
          // print a header to the data file.

          if (NERRS == 0) {
            NOUNIT.println(
              '\n DGV -- Real Generalized eigenvalue problem driver',
            );

            // Matrix types

            NOUNIT.println(' Matrix types (see DDRGEV3 for details): ');
            NOUNIT.println(
              ' Special Matrices:${' ' * 23}(J\'=transposed Jordan block)\n   1=(0,0)  2=(I,0)  3=(0,I)  4=(I,I)  5=(J\',J\')  6=(diag(J\',I), diag(I,J\'))\n Diagonal Matrices:  ( D=diag(0,1,2,...) )\n   7=(D,I)   9=(large*D, small*I)  11=(large*I, small*D)  13=(large*D, large*I)\n   8=(I,D)  10=(small*D, large*I)  12=(small*I, large*D)  14=(small*D, small*I)\n  15=(D, reversed D),',
            );
            NOUNIT.println(
              ' Matrices Rotated by Random Orthogonal Matrices U, V:\n  16=Transposed Jordan Blocks             19=geometric alpha, beta=0,1\n  17=arithm. alpha&beta                   20=arithmetic alpha, beta=0,1\n  18=clustered alpha, beta=0,1            21=random alpha, beta=0,1\n Large & Small Matrices:\n  22=(large, small)   23=(small,large)    24=(small,small)    25=(large,large)\n  26=random O(1) matrices.',
            );

            // Tests performed

            NOUNIT.println(
              '\n Tests performed:    \n 1 = max | ( b A - a B )\'*l | / const.,\n 2 = | |VR(i)| - 1 | / ulp,\n 3 = max | ( b A - a B )*r | / const.\n 4 = | |VL(i)| - 1 | / ulp,\n 5 = 0 if W same no matter if r or l computed,\n 6 = 0 if l same no matter if l computed,\n 7 = 0 if r same no matter if r computed,/n ',
            );
          }
          NERRS = NERRS + 1;
          if (RESULT[JR] < 10000.0) {
            NOUNIT.println(
              ' Matrix order=${N.i5}, type=${JTYPE.i2}, seed=${IOLDSD.i4(4, ',')} result ${JR.i2} is${RESULT[JR].f8_2}',
            );
          } else {
            NOUNIT.println(
              ' Matrix order=${N.i5}, type=${JTYPE.i2}, seed=${IOLDSD.i4(4, ',')} result ${JR.i2} is${(RESULT[JR] * 10).d10_3}',
            );
          }
        }
      }
    }
  }

  // Summary
  alasvm('DGV', NOUNIT, NERRS, NTESTT, 0);

  WORK[1] = MAXWRK.toDouble();
}

void _print9999(
  final Nout NOUNIT,
  final String s,
  final int info,
  final int n,
  final int jtype,
  final Array<int> iseed,
) {
  NOUNIT.println(
    ' DDRGEV3: $s returned INFO=${info.i6}.\n${' ' * 3}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i4(4, ',')})',
  );
}

void _print9998(
  final Nout NOUNIT,
  final String side,
  final String s,
  final double error,
  final int n,
  final int jtype,
  final Array<int> iseed,
) {
  NOUNIT.println(
    ' DDRGEV3: $side Eigenvectors from $s incorrectly normalized.\n Bits of error=${error.g10_3},${' ' * 3}N=${n.i4}, JTYPE=${jtype.i3}, ISEED=(${iseed.i4(4, ',')})',
  );
}
