// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dggev3.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlarfg.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dorm2r.dart';
import 'package:dart_lapack/src/intrinsics/sign.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/range.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:test/test.dart';

import '../matgen/dlarnd.dart';
import '../test_driver.dart';
import 'alasvm.dart';
import 'dget52.dart';
import 'dlatm4.dart';
import 'xlaenv.dart';

void ddrgev3(
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
  final Matrix<double> QE_,
  final int LDQE,
  final Array<double> ALPHAR_,
  final Array<double> ALPHAI_,
  final Array<double> BETA_,
  final Array<double> ALPHR1_,
  final Array<double> ALPHI1_,
  final Array<double> BETA1_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<double> RESULT_,
  final Box<int> INFO,
  final TestDriver test,
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
  final QE = QE_.having(ld: LDQE);
  final ALPHAR = ALPHAR_.having();
  final ALPHAI = ALPHAI_.having();
  final BETA = BETA_.having();
  final ALPHR1 = ALPHR1_.having();
  final ALPHI1 = ALPHI1_.having();
  final BETA1 = BETA1_.having();
  final WORK = WORK_.having();
  final RESULT = RESULT_.having();
  const ZERO = 0.0, ONE = 1.0;
  const MAXTYP = 27;
  const KCLASS = [
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    3, 4 //
  ];
  const KZ1 = [0, 1, 2, 1, 3, 3];
  const KZ2 = [0, 0, 1, 2, 1, 1];
  const KADD = [0, 0, 0, 0, 3, 2];
  const KATYPE = [
    0, 1, 0, 1, 2, 3, 4, 1, 4, 4, 1, 1, 4, 4, //
    4, 2, 4, 5, 8, 7, 9, 4, 4, 4, 4, 0, 0 //
  ];
  const KBTYPE = [
    0, 0, 1, 1, 2, -3, 1, 4, 1, 1, 4, 4, 1, 1, -4, 2, -4, 8, 8, 8, 8, 8, 8, 8,
    8, 0, 0 //
  ];
  const KAZERO = [
    1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 1, 2, 2, 3, 1, 3, 5, 5, 5, 5, 3, 3, 3, 3,
    1, 1 //
  ];
  const KBZERO = [
    1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 4, 1, 4, 6, 6, 6, 6, 4, 4, 4, 4,
    1, 1 //
  ];
  const KAMAGN = [
    1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 2, 3, 2, 3, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 2,
    1, 3 //
  ];
  const KBMAGN = [
    1, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2, 2, 3, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2,
    1, 3 //
  ];
  const KTRIAN = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1 //
  ];
  const IASIGN = [
    0, 0, 0, 0, 0, 0, 2, 0, 2, 2, 0, 0, 2, 2, 2, 0, 2, 0, 0, 0, 2, 2, 2, 2, 2,
    0, 0 //
  ];
  const IBSIGN = [
    0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 2, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0 //
  ];

  // Check for errors

  INFO.value = 0;
  var MAXWRK = 0;
  {
    var BADNN = false;
    var NMAX = 1;
    for (var J = 1; J <= NSIZES; J++) {
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

    var MINWRK = 1;
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
  }

  // Quick return if possible

  if (NSIZES == 0 || NTYPES == 0) return;

  final ULP = dlamch('Epsilon') * dlamch('Base');
  final SAFMIN = dlamch('Safe minimum') / ULP;
  final SAFMAX = ONE / SAFMIN;
  final ULPINV = ONE / ULP;

  // Loop over sizes, types

  var NTESTT = 0;
  var NERRS = 0;

  for (final JSIZE in 1.through(NSIZES)) {
    final N = NN[JSIZE];
    final N1 = max(1, N);

    final RMAGN = Array.fromList([
      ZERO,
      ONE,
      SAFMAX * ULP / N1,
      SAFMIN * ULPINV * N1,
    ], offset: zeroIndexedArrayOffset);

    final MTYPES = NSIZES != 1 ? min(MAXTYP, NTYPES) : min(MAXTYP + 1, NTYPES);

    for (final JTYPE in 1.through(MTYPES)) {
      final skip = !DOTYPE[JTYPE];
      test('DDRGEV3 (SIZE = $N, TYPE = $JTYPE)', () {
        // Save ISEED in case of an error.
        final IOLDSD = ISEED.copy();
        final IERR = Box(0);

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
          if (KCLASS[JTYPE - 1] < 3) {
            // Generate A (w/o rotation)

            int IN;
            if (KATYPE[JTYPE - 1].abs() == 3) {
              IN = 2 * ((N - 1) ~/ 2) + 1;
              if (IN != N) dlaset('Full', N, N, ZERO, ZERO, A, LDA);
            } else {
              IN = N;
            }
            dlatm4(
                KATYPE[JTYPE - 1],
                IN,
                KZ1[KAZERO[JTYPE - 1] - 1],
                KZ2[KAZERO[JTYPE - 1] - 1],
                IASIGN[JTYPE - 1],
                RMAGN[KAMAGN[JTYPE - 1]],
                ULP,
                RMAGN[KTRIAN[JTYPE - 1] * KAMAGN[JTYPE - 1]],
                2,
                ISEED,
                A,
                LDA);
            var IADD = KADD[KAZERO[JTYPE - 1] - 1];
            if (IADD > 0 && IADD <= N) A[IADD][IADD] = ONE;

            // Generate B (w/o rotation)

            if (KBTYPE[JTYPE - 1].abs() == 3) {
              IN = 2 * ((N - 1) ~/ 2) + 1;
              if (IN != N) dlaset('Full', N, N, ZERO, ZERO, B, LDA);
            } else {
              IN = N;
            }
            dlatm4(
                KBTYPE[JTYPE - 1],
                IN,
                KZ1[KBZERO[JTYPE - 1] - 1],
                KZ2[KBZERO[JTYPE - 1] - 1],
                IBSIGN[JTYPE - 1],
                RMAGN[KBMAGN[JTYPE - 1]],
                ONE,
                RMAGN[KTRIAN[JTYPE - 1] * KBMAGN[JTYPE - 1]],
                2,
                ISEED,
                B,
                LDA);
            IADD = KADD[KBZERO[JTYPE - 1] - 1];
            if (IADD != 0 && IADD <= N) B[IADD][IADD] = ONE;

            if (KCLASS[JTYPE - 1] == 2 && N > 0) {
              // Include rotations

              // Generate Q, Z as Householder transformations times
              // a diagonal matrix.

              for (var JC = 1; JC <= N - 1; JC++) {
                for (var JR = JC; JR <= N; JR++) {
                  Q[JR][JC] = dlarnd(3, ISEED);
                  Z[JR][JC] = dlarnd(3, ISEED);
                }
                dlarfg(N + 1 - JC, Q.box(JC, JC), Q(JC + 1, JC).asArray(), 1,
                    WORK.box(JC));
                WORK[2 * N + JC] = sign(ONE, Q[JC][JC]);
                Q[JC][JC] = ONE;
                dlarfg(N + 1 - JC, Z.box(JC, JC), Z(JC + 1, JC).asArray(), 1,
                    WORK.box(N + JC));
                WORK[3 * N + JC] = sign(ONE, Z[JC][JC]);
                Z[JC][JC] = ONE;
              }
              Q[N][N] = ONE;
              WORK[N] = ZERO;
              WORK[3 * N] = sign(ONE, dlarnd(2, ISEED));
              Z[N][N] = ONE;
              WORK[2 * N] = ZERO;
              WORK[4 * N] = sign(ONE, dlarnd(2, ISEED));

              // Apply the diagonal matrices

              for (var JC = 1; JC <= N; JC++) {
                for (var JR = 1; JR <= N; JR++) {
                  A[JR][JC] = WORK[2 * N + JR] * WORK[3 * N + JC] * A[JR][JC];
                  B[JR][JC] = WORK[2 * N + JR] * WORK[3 * N + JC] * B[JR][JC];
                }
              }
              dorm2r('L', 'N', N, N, N - 1, Q, LDQ, WORK, A, LDA,
                  WORK(2 * N + 1), IERR);
              if (IERR.value == 0) {
                dorm2r('R', 'T', N, N, N - 1, Z, LDQ, WORK(N + 1), A, LDA,
                    WORK(2 * N + 1), IERR);
                if (IERR.value == 0) {
                  dorm2r('L', 'N', N, N, N - 1, Q, LDQ, WORK, B, LDA,
                      WORK(2 * N + 1), IERR);
                  if (IERR.value == 0) {
                    dorm2r('R', 'T', N, N, N - 1, Z, LDQ, WORK(N + 1), B, LDA,
                        WORK(2 * N + 1), IERR);
                  }
                }
              }
            }
          } else if (KCLASS[JTYPE - 1] == 3) {
            // Random matrices

            for (var JC = 1; JC <= N; JC++) {
              for (var JR = 1; JR <= N; JR++) {
                A[JR][JC] = RMAGN[KAMAGN[JTYPE - 1]] * dlarnd(2, ISEED);
                B[JR][JC] = RMAGN[KBMAGN[JTYPE - 1]] * dlarnd(2, ISEED);
              }
            }
          } else {
            // Random upper Hessenberg pencil with singular B

            for (var JC = 1; JC <= N; JC++) {
              for (var JR = 1; JR <= min(JC + 1, N); JR++) {
                A[JR][JC] = RMAGN[KAMAGN[JTYPE - 1]] * dlarnd(2, ISEED);
              }
              for (var JR = JC + 2; JR <= N; JR++) {
                A[JR][JC] = ZERO;
              }
            }
            for (var JC = 1; JC <= N; JC++) {
              for (var JR = 1; JR <= JC; JR++) {
                B[JR][JC] = RMAGN[KAMAGN[JTYPE - 1]] * dlarnd(2, ISEED);
              }
              for (var JR = JC + 1; JR <= N; JR++) {
                B[JR][JC] = ZERO;
              }
            }
            for (var JC = 1; JC <= N; JC += 4) {
              B[JC][JC] = ZERO;
            }
          }

          if (IERR.value != 0) {
            _print9999(NOUNIT, 'Generator', IERR.value, N, JTYPE, IOLDSD);
            INFO.value = (IERR.value).abs();
            return;
          }
        }

        for (var I = 1; I <= 7; I++) {
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
          dggev3('V', 'V', N, S, LDA, T, LDA, ALPHAR, ALPHAI, BETA, Q, LDQ, Z,
              LDQ, WORK, LWORK, IERR);
          if (IERR.value != 0 && IERR.value != N + 1) {
            RESULT[1] = ULPINV;
            _print9999(NOUNIT, 'DGGEV31', IERR.value, N, JTYPE, IOLDSD);
            INFO.value = (IERR.value).abs();
            break;
          }

          // Do the tests (1) and (2)

          dget52(true, N, A, LDA, B, LDA, Q, LDQ, ALPHAR, ALPHAI, BETA, WORK,
              RESULT(1));
          if (RESULT[2] > THRESH) {
            _print9998(NOUNIT, 'Left', 'DGGEV31', RESULT[2], N, JTYPE, IOLDSD);
          }

          // Do the tests (3) and (4)

          dget52(false, N, A, LDA, B, LDA, Z, LDQ, ALPHAR, ALPHAI, BETA, WORK,
              RESULT(3));
          if (RESULT[4] > THRESH) {
            _print9998(NOUNIT, 'Right', 'DGGEV31', RESULT[4], N, JTYPE, IOLDSD);
          }

          // Do the test (5)

          dlacpy(' ', N, N, A, LDA, S, LDA);
          dlacpy(' ', N, N, B, LDA, T, LDA);
          dggev3('N', 'N', N, S, LDA, T, LDA, ALPHR1, ALPHI1, BETA1, Q, LDQ, Z,
              LDQ, WORK, LWORK, IERR);
          if (IERR.value != 0 && IERR.value != N + 1) {
            RESULT[1] = ULPINV;
            _print9999(NOUNIT, 'DGGEV32', IERR.value, N, JTYPE, IOLDSD);
            INFO.value = (IERR.value).abs();
            break;
          }

          for (var J = 1; J <= N; J++) {
            if (ALPHAR[J] != ALPHR1[J] ||
                ALPHAI[J] != ALPHI1[J] ||
                BETA[J] != BETA1[J]) RESULT[5] = ULPINV;
          }

          // Do the test (6): Compute eigenvalues and left eigenvectors,
          // and test them

          dlacpy(' ', N, N, A, LDA, S, LDA);
          dlacpy(' ', N, N, B, LDA, T, LDA);
          dggev3('V', 'N', N, S, LDA, T, LDA, ALPHR1, ALPHI1, BETA1, QE, LDQE,
              Z, LDQ, WORK, LWORK, IERR);
          if (IERR.value != 0 && IERR.value != N + 1) {
            RESULT[1] = ULPINV;
            _print9999(NOUNIT, 'DGGEV33', IERR.value, N, JTYPE, IOLDSD);
            INFO.value = (IERR.value).abs();
            break;
          }

          for (var J = 1; J <= N; J++) {
            if (ALPHAR[J] != ALPHR1[J] ||
                ALPHAI[J] != ALPHI1[J] ||
                BETA[J] != BETA1[J]) RESULT[6] = ULPINV;
          }

          for (var J = 1; J <= N; J++) {
            for (var JC = 1; JC <= N; JC++) {
              if (Q[J][JC] != QE[J][JC]) RESULT[6] = ULPINV;
            }
          }

          // DO the test (7): Compute eigenvalues and right eigenvectors,
          // and test them

          dlacpy(' ', N, N, A, LDA, S, LDA);
          dlacpy(' ', N, N, B, LDA, T, LDA);
          dggev3('N', 'V', N, S, LDA, T, LDA, ALPHR1, ALPHI1, BETA1, Q, LDQ, QE,
              LDQE, WORK, LWORK, IERR);
          if (IERR.value != 0 && IERR.value != N + 1) {
            RESULT[1] = ULPINV;
            _print9999(NOUNIT, 'DGGEV34', IERR.value, N, JTYPE, IOLDSD);
            INFO.value = (IERR.value).abs();
            break;
          }

          for (var J = 1; J <= N; J++) {
            if (ALPHAR[J] != ALPHR1[J] ||
                ALPHAI[J] != ALPHI1[J] ||
                BETA[J] != BETA1[J]) RESULT[7] = ULPINV;
          }

          for (var J = 1; J <= N; J++) {
            for (var JC = 1; JC <= N; JC++) {
              if (Z[J][JC] != QE[J][JC]) RESULT[7] = ULPINV;
            }
          }

          // End of Loop -- Check for RESULT[j] > THRESH
          break;
        }

        NTESTT += 7;

        // Print out tests which fail.

        for (var JR = 1; JR <= 7; JR++) {
          final reason = RESULT[JR] < 10000.0
              ? ' Matrix order=${N.i5}, type=${JTYPE.i2}, seed=${IOLDSD.i4(4, ',')} result ${JR.i2} is${RESULT[JR].f8_2}'
              : ' Matrix order=${N.i5}, type=${JTYPE.i2}, seed=${IOLDSD.i4(4, ',')} result ${JR.i2} is${(RESULT[JR] * 10).d10_3}';
          test.expect(RESULT[JR], lessThan(THRESH), reason: reason);
          if (RESULT[JR] >= THRESH) {
            // If this is the first test to fail,
            // print a header to the data file.

            if (NERRS == 0) {
              NOUNIT.println(
                  '\n DGV -- Real Generalized eigenvalue problem driver');

              // Matrix types

              NOUNIT.println(' Matrix types (see DDRGEV3 for details): ');
              NOUNIT.println(
                  ' Special Matrices:${' ' * 23}(J\'=transposed Jordan block)\n   1=(0,0)  2=(I,0)  3=(0,I)  4=(I,I)  5=(J\',J\')  6=(diag(J\',I), diag(I,J\'))\n Diagonal Matrices:  ( D=diag(0,1,2,...) )\n   7=(D,I)   9=(large*D, small*I)  11=(large*I, small*D)  13=(large*D, large*I)\n   8=(I,D)  10=(small*D, large*I)  12=(small*I, large*D)  14=(small*D, small*I)\n  15=(D, reversed D),');
              NOUNIT.println(
                  ' Matrices Rotated by Random Orthogonal Matrices U, V:\n  16=Transposed Jordan Blocks             19=geometric alpha, beta=0,1\n  17=arithm. alpha&beta                   20=arithmetic alpha, beta=0,1\n  18=clustered alpha, beta=0,1            21=random alpha, beta=0,1\n Large & Small Matrices:\n  22=(large, small)   23=(small,large)    24=(small,small)    25=(large,large)\n  26=random O(1) matrices.');

              // Tests performed

              NOUNIT.println(
                  '\n Tests performed:    \n 1 = max | ( b A - a B )\'*l | / const.,\n 2 = | |VR(i)| - 1 | / ulp,\n 3 = max | ( b A - a B )*r | / const.\n 4 = | |VL(i)| - 1 | / ulp,\n 5 = 0 if W same no matter if r or l computed,\n 6 = 0 if l same no matter if l computed,\n 7 = 0 if r same no matter if r computed,/n ');
            }
            NERRS++;
            NOUNIT.println(reason);
          }
        }
      }, skip: skip);
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
      ' DDRGEV3: $s returned INFO=${info.i6}.\n   N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i4(4, ',')})');
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
      ' DDRGEV3: $side Eigenvectors from $s incorrectly normalized.\n Bits of error=${error.g10_3},   N=${n.i4}, JTYPE=${jtype.i3}, ISEED=(${iseed.i4(4, ',')})');
}
