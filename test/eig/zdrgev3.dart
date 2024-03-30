import 'dart:math';

import 'package:collection/collection.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zggev3.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlarfg.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zunm2r.dart';

import '../matgen/zlarnd.dart';
import 'alasvm.dart';
import 'xlaenv.dart';
import 'zget52.dart';
import 'zlatm4.dart';

void zdrgev3(
  final int NSIZES,
  final Array<int> NN_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final Array<int> ISEED_,
  final double THRESH,
  final Nout NOUNIT,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final Matrix<Complex> S_,
  final Matrix<Complex> T_,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Matrix<Complex> Z_,
  final Matrix<Complex> QE_,
  final int LDQE,
  final Array<Complex> ALPHA_,
  final Array<Complex> BETA_,
  final Array<Complex> ALPHA1_,
  final Array<Complex> BETA1_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<double> RESULT_,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final NN = NN_.having();
  final DOTYPE = DOTYPE_.having();
  final ISEED = ISEED_.having(length: 4);
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDA);
  final S = S_.having(ld: LDA);
  final T = T_.having(ld: LDA);
  final Q = Q_.having(ld: LDQ);
  final QE = QE_.having(ld: LDQE);
  final Z = Z_.having(ld: LDQ);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final ALPHA = ALPHA_.having();
  final BETA = BETA_.having();
  final ALPHA1 = ALPHA1_.having();
  final BETA1 = BETA1_.having();
  final RESULT = RESULT_.having();
  const ZERO = 0.0, ONE = 1.0;
  const MAXTYP = 26;
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
      N = 0,
      N1,
      NB,
      NERRS,
      NMAX,
      NTESTT;
  double SAFMAX, SAFMIN, ULP, ULPINV;
  Complex CTEMP;
  final IOLDSD = Array<int>(4);
  final RMAGN = Array<double>(4, offset: zeroIndexedArrayOffset);
  final KCLASS = Array.fromList([
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, //
    1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3
  ]);
  final KZ1 = Array.fromList([0, 1, 2, 1, 3, 3]);
  final KZ2 = Array.fromList([0, 0, 1, 2, 1, 1]);
  final KADD = Array.fromList([0, 0, 0, 0, 3, 2]);
  final KATYPE = Array.fromList([
    0, 1, 0, 1, 2, 3, 4, 1, 4, 4, 1, 1, 4, //
    4, 4, 2, 4, 5, 8, 7, 9, 4, 4, 4, 4, 0
  ]);
  final KBTYPE = Array.fromList([
    0, 0, 1, 1, 2, -3, 1, 4, 1, 1, 4, 4, 1, //
    1, -4, 2, -4, 8, 8, 8, 8, 8, 8, 8, 8, 0
  ]);
  final KAZERO = Array.fromList([
    1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 1, 2, //
    2, 3, 1, 3, 5, 5, 5, 5, 3, 3, 3, 3, 1
  ]);
  final KBZERO = Array.fromList([
    1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 1, //
    1, 4, 1, 4, 6, 6, 6, 6, 4, 4, 4, 4, 1
  ]);
  final KAMAGN = Array.fromList([
    1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 2, 3, 2, //
    3, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 2, 1
  ]);
  final KBMAGN = Array.fromList([
    1, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2, 2, //
    3, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2, 1
  ]);
  final KTRIAN = Array.fromList([
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //
    0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
  ]);
  final LASIGN = Array.fromList([
    false, false, false, false, false, false, true, false, true, true, //
    false, false, true, true, true, false, true, false, false, false, //
    true, true, true, true, true, false
  ]);
  final LBSIGN = Array.fromList([
    false, false, false, false, false, false, false, true, false, false, //
    true, true, false, false, true, false, true, false, false, false, //
    false, false, false, false, false, false
  ]);
  final IERR = Box(0);

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
  //  (Note: Comments in the code beginning "Workspace:" describe the
  //   minimal amount of workspace needed at that point in the code,
  //   as well as the preferred amount for good performance.
  //   NB refers to the optimal block size for the immediately
  //   following subroutine, as returned by ILAENV.

  MINWRK = 1;
  if (INFO.value == 0 && LWORK >= 1) {
    MINWRK = NMAX * (NMAX + 1);
    NB = [
      1,
      ilaenv(1, 'ZGEQRF', ' ', NMAX, NMAX, -1, -1),
      ilaenv(1, 'ZUNMQR', 'LC', NMAX, NMAX, NMAX, -1),
      ilaenv(1, 'ZUNGQR', ' ', NMAX, NMAX, NMAX, -1)
    ].max;
    MAXWRK = [2 * NMAX, NMAX * (NB + 1), NMAX * (NMAX + 1)].max;
    WORK[1] = MAXWRK.toComplex();
  }

  if (LWORK < MINWRK) INFO.value = -23;

  if (INFO.value != 0) {
    xerbla('ZDRGEV3', -INFO.value);
    return;
  }

  // Quick return if possible

  if (NSIZES == 0 || NTYPES == 0) return;

  ULP = dlamch('Precision');
  SAFMIN = dlamch('Safe minimum');
  SAFMIN /= ULP;
  SAFMAX = ONE / SAFMIN;
  ULPINV = ONE / ULP;

  // The values RMAGN(2:3) depend on N, see below.

  RMAGN[0] = ZERO;
  RMAGN[1] = ONE;

  // Loop over sizes, types

  NTESTT = 0;
  NERRS = 0;

  for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) {
    N = NN[JSIZE];
    N1 = max(1, N);
    RMAGN[2] = SAFMAX * ULP / N1;
    RMAGN[3] = SAFMIN * ULPINV * N1;

    if (NSIZES != 1) {
      MTYPES = min(MAXTYP, NTYPES);
    } else {
      MTYPES = min(MAXTYP + 1, NTYPES);
    }

    for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
      if (!DOTYPE[JTYPE]) continue;

      // Save ISEED in case of an error.

      for (J = 1; J <= 4; J++) {
        IOLDSD[J] = ISEED[J];
      }

      // Generate test matrices A and B

      // Description of control parameters:

      // KZLASS: =1 means w/o rotation, =2 means w/ rotation,
      //         =3 means random.
      // KATYPE: the "type" to be passed to ZLATM4 for computing A.
      // KAZERO: the pattern of zeros on the diagonal for A:
      //         =1: ( xxx ), =2: (0, xxx ) =3: ( 0, 0, xxx, 0 ),
      //         =4: ( 0, xxx, 0, 0 ), =5: ( 0, 0, 1, xxx, 0 ),
      //         =6: ( 0, 1, 0, xxx, 0 ).  (xxx means a string of
      //         non-zero entries.)
      // KAMAGN: the magnitude of the matrix: =0: zero, =1: O(1),
      //         =2: large, =3: small.
      // LASIGN: true if the diagonal elements of A are to be
      //         multiplied by a random magnitude 1 number.
      // KBTYPE, KBZERO, KBMAGN, LBSIGN: the same, but for B.
      // KTRIAN: =0: don't fill in the upper triangle, =1: do.
      // KZ1, KZ2, KADD: used to implement KAZERO and KBZERO.
      // RMAGN: used to implement KAMAGN and KBMAGN.

      if (MTYPES <= MAXTYP) {
        IERR.value = 0;
        if (KCLASS[JTYPE] < 3) {
          // Generate A (w/o rotation)

          if (KATYPE[JTYPE].abs() == 3) {
            IN = 2 * ((N - 1) ~/ 2) + 1;
            if (IN != N) {
              zlaset('Full', N, N, Complex.zero, Complex.zero, A, LDA);
            }
          } else {
            IN = N;
          }
          zlatm4(
              KATYPE[JTYPE],
              IN,
              KZ1[KAZERO[JTYPE]],
              KZ2[KAZERO[JTYPE]],
              LASIGN[JTYPE],
              RMAGN[KAMAGN[JTYPE]],
              ULP,
              RMAGN[KTRIAN[JTYPE] * KAMAGN[JTYPE]],
              2,
              ISEED,
              A,
              LDA);
          IADD = KADD[KAZERO[JTYPE]];
          if (IADD > 0 && IADD <= N) {
            A[IADD][IADD] = RMAGN[KAMAGN[JTYPE]].toComplex();
          }

          // Generate B (w/o rotation)

          if (KBTYPE[JTYPE].abs() == 3) {
            IN = 2 * ((N - 1) ~/ 2) + 1;
            if (IN != N) {
              zlaset('Full', N, N, Complex.zero, Complex.zero, B, LDA);
            }
          } else {
            IN = N;
          }
          zlatm4(
              KBTYPE[JTYPE],
              IN,
              KZ1[KBZERO[JTYPE]],
              KZ2[KBZERO[JTYPE]],
              LBSIGN[JTYPE],
              RMAGN[KBMAGN[JTYPE]],
              ONE,
              RMAGN[KTRIAN[JTYPE] * KBMAGN[JTYPE]],
              2,
              ISEED,
              B,
              LDA);
          IADD = KADD[KBZERO[JTYPE]];
          if (IADD != 0 && IADD <= N) {
            B[IADD][IADD] = RMAGN[KBMAGN[JTYPE]].toComplex();
          }

          if (KCLASS[JTYPE] == 2 && N > 0) {
            // Include rotations

            // Generate Q, Z as Householder transformations times
            // a diagonal matrix.

            for (JC = 1; JC <= N - 1; JC++) {
              for (JR = JC; JR <= N; JR++) {
                Q[JR][JC] = zlarnd(3, ISEED);
                Z[JR][JC] = zlarnd(3, ISEED);
              }
              zlarfg(
                  N + 1 - JC, Q(JC, JC), Q(JC + 1, JC).asArray(), 1, WORK(JC));
              WORK[2 * N + JC] = sign(ONE, Q[JC][JC].real).toComplex();
              Q[JC][JC] = Complex.one;
              zlarfg(N + 1 - JC, Z(JC, JC), Z(JC + 1, JC).asArray(), 1,
                  WORK(N + JC));
              WORK[3 * N + JC] = sign(ONE, Z[JC][JC].real).toComplex();
              Z[JC][JC] = Complex.one;
            }
            CTEMP = zlarnd(3, ISEED);
            Q[N][N] = Complex.one;
            WORK[N] = Complex.zero;
            WORK[3 * N] = CTEMP / CTEMP.abs().toComplex();
            CTEMP = zlarnd(3, ISEED);
            Z[N][N] = Complex.one;
            WORK[2 * N] = Complex.zero;
            WORK[4 * N] = CTEMP / CTEMP.abs().toComplex();

            // Apply the diagonal matrices

            for (JC = 1; JC <= N; JC++) {
              for (JR = 1; JR <= N; JR++) {
                A[JR][JC] =
                    WORK[2 * N + JR] * WORK[3 * N + JC].conjugate() * A[JR][JC];
                B[JR][JC] =
                    WORK[2 * N + JR] * WORK[3 * N + JC].conjugate() * B[JR][JC];
              }
            }
            zunm2r('L', 'N', N, N, N - 1, Q, LDQ, WORK, A, LDA, WORK(2 * N + 1),
                IERR);
            if (IERR.value == 0) {
              zunm2r('R', 'C', N, N, N - 1, Z, LDQ, WORK(N + 1), A, LDA,
                  WORK(2 * N + 1), IERR);
              if (IERR.value == 0) {
                zunm2r('L', 'N', N, N, N - 1, Q, LDQ, WORK, B, LDA,
                    WORK(2 * N + 1), IERR);
                if (IERR.value == 0) {
                  zunm2r('R', 'C', N, N, N - 1, Z, LDQ, WORK(N + 1), B, LDA,
                      WORK(2 * N + 1), IERR);
                }
              }
            }
          }
        } else {
          // Random matrices

          for (JC = 1; JC <= N; JC++) {
            for (JR = 1; JR <= N; JR++) {
              A[JR][JC] = RMAGN[KAMAGN[JTYPE]].toComplex() * zlarnd(4, ISEED);
              B[JR][JC] = RMAGN[KBMAGN[JTYPE]].toComplex() * zlarnd(4, ISEED);
            }
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

      // Call XLAENV to set the parameters used in ZLAQZ0

      xlaenv(12, 10);
      xlaenv(13, 12);
      xlaenv(14, 13);
      xlaenv(15, 2);
      xlaenv(17, 10);

      while (true) {
        // Call ZGGEV3 to compute eigenvalues and eigenvectors.

        zlacpy(' ', N, N, A, LDA, S, LDA);
        zlacpy(' ', N, N, B, LDA, T, LDA);
        zggev3('V', 'V', N, S, LDA, T, LDA, ALPHA, BETA, Q, LDQ, Z, LDQ, WORK,
            LWORK, RWORK, IERR);
        if (IERR.value != 0 && IERR.value != N + 1) {
          RESULT[1] = ULPINV;
          _print9999(NOUNIT, 'ZGGEV31', IERR.value, N, JTYPE, IOLDSD);
          INFO.value = (IERR.value).abs();
          break;
        }

        // Do the tests (1) and (2)

        zget52(true, N, A, LDA, B, LDA, Q, LDQ, ALPHA, BETA, WORK, RWORK,
            RESULT(1));
        if (RESULT[2] > THRESH) {
          _print9998(NOUNIT, 'Left', 'ZGGEV31', RESULT[2], N, JTYPE, IOLDSD);
        }

        // Do the tests (3) and (4)

        zget52(false, N, A, LDA, B, LDA, Z, LDQ, ALPHA, BETA, WORK, RWORK,
            RESULT(3));
        if (RESULT[4] > THRESH) {
          _print9998(NOUNIT, 'Right', 'ZGGEV31', RESULT[4], N, JTYPE, IOLDSD);
        }

        // Do test (5)

        zlacpy(' ', N, N, A, LDA, S, LDA);
        zlacpy(' ', N, N, B, LDA, T, LDA);
        zggev3('N', 'N', N, S, LDA, T, LDA, ALPHA1, BETA1, Q, LDQ, Z, LDQ, WORK,
            LWORK, RWORK, IERR);
        if (IERR.value != 0 && IERR.value != N + 1) {
          RESULT[1] = ULPINV;
          _print9999(NOUNIT, 'ZGGEV32', IERR.value, N, JTYPE, IOLDSD);
          INFO.value = (IERR.value).abs();
          break;
        }

        for (J = 1; J <= N; J++) {
          if (ALPHA[J] != ALPHA1[J] || BETA[J] != BETA1[J]) RESULT[5] = ULPINV;
        }

        // Do test (6): Compute eigenvalues and left eigenvectors,
        // and test them

        zlacpy(' ', N, N, A, LDA, S, LDA);
        zlacpy(' ', N, N, B, LDA, T, LDA);
        zggev3('V', 'N', N, S, LDA, T, LDA, ALPHA1, BETA1, QE, LDQE, Z, LDQ,
            WORK, LWORK, RWORK, IERR);
        if (IERR.value != 0 && IERR.value != N + 1) {
          RESULT[1] = ULPINV;
          _print9999(NOUNIT, 'ZGGEV33', IERR.value, N, JTYPE, IOLDSD);
          INFO.value = (IERR.value).abs();
          break;
        }

        for (J = 1; J <= N; J++) {
          if (ALPHA[J] != ALPHA1[J] || BETA[J] != BETA1[J]) RESULT[6] = ULPINV;
        }

        for (J = 1; J <= N; J++) {
          for (JC = 1; JC <= N; JC++) {
            if (Q(J, JC) != QE(J, JC)) RESULT[6] = ULPINV;
          }
        }

        // Do test (7): Compute eigenvalues and right eigenvectors,
        // and test them

        zlacpy(' ', N, N, A, LDA, S, LDA);
        zlacpy(' ', N, N, B, LDA, T, LDA);
        zggev3('N', 'V', N, S, LDA, T, LDA, ALPHA1, BETA1, Q, LDQ, QE, LDQE,
            WORK, LWORK, RWORK, IERR);
        if (IERR.value != 0 && IERR.value != N + 1) {
          RESULT[1] = ULPINV;
          _print9999(NOUNIT, 'ZGGEV34', IERR.value, N, JTYPE, IOLDSD);
          INFO.value = (IERR.value).abs();
          break;
        }

        for (J = 1; J <= N; J++) {
          if (ALPHA[J] != ALPHA1[J] || BETA[J] != BETA1[J]) RESULT[7] = ULPINV;
        }

        for (J = 1; J <= N; J++) {
          for (JC = 1; JC <= N; JC++) {
            if (Z(J, JC) != QE(J, JC)) RESULT[7] = ULPINV;
          }
        }

        break;
      }

      // End of Loop -- Check for RESULT(j) > THRESH

      NTESTT += 7;

      // Print out tests which fail.

      for (JR = 1; JR <= 7; JR++) {
        if (RESULT[JR] >= THRESH) {
          // If this is the first test to fail,
          // print a header to the data file.

          if (NERRS == 0) {
            NOUNIT.println(
                '\n ZGV -- Complex Generalized eigenvalue problem driver');

            // Matrix types

            NOUNIT.println(' Matrix types (see ZDRGEV3 for details): ');
            NOUNIT.println(' Special Matrices:${' ' * 23}(J'
                '=transposed Jordan block)\n   1=(0,0)  2=(I,0)  3=(0,I)  4=(I,I)  5=(J'
                ',J'
                ')  6=(diag(J'
                ',I), diag(I,J'
                '))\n Diagonal Matrices:  ( D=diag(0,1,2,...) )\n   7=(D,I)   9=(large*D, small*I)  11=(large*I, small*D)  13=(large*D, large*I)\n   8=(I,D)  10=(small*D, large*I)  12=(small*I, large*D)  14=(small*D, small*I)\n  15=(D, reversed D)');
            NOUNIT.println(
                ' Matrices Rotated by Random Orthogonal Matrices U, V:\n  16=Transposed Jordan Blocks             19=geometric alpha, beta=0,1\n  17=arithm. alpha&beta                   20=arithmetic alpha, beta=0,1\n  18=clustered alpha, beta=0,1            21=random alpha, beta=0,1\n Large & Small Matrices:\n  22=(large, small)   23=(small,large)    24=(small,small)    25=(large,large)\n  26=random O(1) matrices.');

            // Tests performed

            NOUNIT.println('\n Tests performed:    \n 1 = max | ( b A - a B )'
                '*l | / const.,\n 2 = | |VR(i)| - 1 | / ulp,\n 3 = max | ( b A - a B )*r | / const.\n 4 = | |VL(i)| - 1 | / ulp,\n 5 = 0 if W same no matter if r or l computed,\n 6 = 0 if l same no matter if l computed,\n 7 = 0 if r same no matter if r computed,/n ');
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

  alasvm('ZGV3', NOUNIT, NERRS, NTESTT, 0);

  WORK[1] = MAXWRK.toComplex();
}

void _print9999(
    Nout nout, String s, int info, int n, int jtype, Array<int> iseed) {
  nout.println(
      ' ZDRGEV3: $s returned INFO=${info.i6}.\n   N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}

void _print9998(Nout nout, String side, String fn, double error, int n,
    int jtype, Array<int> iseed) {
  nout.println(
      ' ZDRGEV3: $side Eigenvectors from $fn incorrectly normalized.\n Bits of error=${error.g10_3},   N=${n.i4}, JTYPE=${jtype.i3}, ISEED=(${iseed.i4(3, ',')})');
}
