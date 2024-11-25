// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dsbtrd.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/xerbla.dart';

import '../matgen/dlatmr.dart';
import '../matgen/dlatms.dart';
import 'dlasum.dart';
import 'dsbt21.dart';

void dchksb(
  final int NSIZES,
  final Array<int> NN_,
  final int NWDTHS,
  final Array<int> KK_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final Array<int> ISEED_,
  final int THRESH,
  final Nout NOUNIT,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> SD_,
  final Array<double> SE_,
  final Matrix<double> U_,
  final int LDU,
  final Array<double> WORK_,
  final int LWORK,
  final Array<double> RESULT_,
  final Box<int> INFO,
) {
  final NN = NN_.having();
  final KK = KK_.having();
  final DOTYPE = DOTYPE_.having();
  final ISEED = ISEED_.having();
  final A = A_.having(ld: LDA);
  final SD = SD_.having();
  final SE = SE_.having();
  final U = U_.having(ld: LDU);
  final WORK = WORK_.having();
  final RESULT = RESULT_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, TEN = 10.0;
  const HALF = ONE / TWO;
  const MAXTYP = 15;
  bool BADNN, BADNNB;
  int I,
      IMODE,
      ITYPE,
      J,
      JC,
      JCOL,
      JR,
      JSIZE,
      JTYPE = 0,
      JWIDTH,
      K = 0,
      KMAX,
      MTYPES,
      N,
      NERRS,
      NMAX,
      NTEST = 0,
      NTESTT;
  double ANINV, ANORM = 0, COND, OVFL, RTOVFL, RTUNFL, TEMP1, ULP, ULPINV, UNFL;
  final IDUMMA = Array<int>(1), IOLDSD = Array<int>(4);
  final IINFO = Box(0);
  final KTYPE = Array.fromList([
    1, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8, //
  ]);
  final KMAGN = Array.fromList([
    1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, //
  ]);
  final KMODE = Array.fromList([
    0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, //
  ]);

  // Check for errors

  NTESTT = 0;
  INFO.value = 0;

  // Important constants

  BADNN = false;
  NMAX = 1;
  for (J = 1; J <= NSIZES; J++) {
    NMAX = max(NMAX, NN[J]);
    if (NN[J] < 0) BADNN = true;
  }

  BADNNB = false;
  KMAX = 0;
  for (J = 1; J <= NSIZES; J++) {
    KMAX = max(KMAX, KK[J]);
    if (KK[J] < 0) BADNNB = true;
  }
  KMAX = min(NMAX - 1, KMAX);

  // Check for errors

  if (NSIZES < 0) {
    INFO.value = -1;
  } else if (BADNN) {
    INFO.value = -2;
  } else if (NWDTHS < 0) {
    INFO.value = -3;
  } else if (BADNNB) {
    INFO.value = -4;
  } else if (NTYPES < 0) {
    INFO.value = -5;
  } else if (LDA < KMAX + 1) {
    INFO.value = -11;
  } else if (LDU < NMAX) {
    INFO.value = -15;
  } else if ((max(LDA, NMAX) + 1) * NMAX > LWORK) {
    INFO.value = -17;
  }

  if (INFO.value != 0) {
    xerbla('DCHKSB', -INFO.value);
    return;
  }

  // Quick return if possible

  if (NSIZES == 0 || NTYPES == 0 || NWDTHS == 0) return;

  // More Important constants

  UNFL = dlamch('Safe minimum');
  OVFL = ONE / UNFL;
  ULP = dlamch('Epsilon') * dlamch('Base');
  ULPINV = ONE / ULP;
  RTUNFL = sqrt(UNFL);
  RTOVFL = sqrt(OVFL);

  // Loop over sizes, types

  NERRS = 0;

  for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) {
    N = NN[JSIZE];
    ANINV = ONE / max(1, N);

    for (JWIDTH = 1; JWIDTH <= NWDTHS; JWIDTH++) {
      K = KK[JWIDTH];
      if (K > N) continue;
      K = max(0, min(N - 1, K));

      if (NSIZES != 1) {
        MTYPES = min(MAXTYP, NTYPES);
      } else {
        MTYPES = min(MAXTYP + 1, NTYPES);
      }

      for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
        if (!DOTYPE[JTYPE]) continue;
        NTEST = 0;

        for (J = 1; J <= 4; J++) {
          IOLDSD[J] = ISEED[J];
        }

        // Compute "A".
        // Store as "Upper"; later, we will copy to other format.
        //
        // Control parameters:
        //
        //     KMAGN  KMODE        KTYPE
        // =1  O(1)   clustered 1  zero
        // =2  large  clustered 2  identity
        // =3  small  exponential  (none)
        // =4         arithmetic   diagonal, (w/ eigenvalues)
        // =5         random log   symmetric, w/ eigenvalues
        // =6         random       (none)
        // =7                      random diagonal
        // =8                      random symmetric
        // =9                      positive definite
        // =10                     diagonally dominant tridiagonal

        if (MTYPES <= MAXTYP) {
          ITYPE = KTYPE[JTYPE];
          IMODE = KMODE[JTYPE];

          // Compute norm

          //  GOTO( 40, 50, 60 )KMAGN[ JTYPE ];
          switch (KMAGN[JTYPE]) {
            case 1:
              ANORM = ONE;
              break;

            case 2:
              ANORM = (RTOVFL * ULP) * ANINV;
              break;

            case 3:
              ANORM = RTUNFL * N * ULPINV;
              break;
          }

          dlaset('Full', LDA, N, ZERO, ZERO, A, LDA);
          IINFO.value = 0;
          if (JTYPE <= 15) {
            COND = ULPINV;
          } else {
            COND = ULPINV * ANINV / TEN;
          }

          // Special Matrices -- Identity & Jordan block

          // Zero

          if (ITYPE == 1) {
            IINFO.value = 0;
          } else if (ITYPE == 2) {
            // Identity

            for (JCOL = 1; JCOL <= N; JCOL++) {
              A[K + 1][JCOL] = ANORM;
            }
          } else if (ITYPE == 4) {
            // Diagonal Matrix, [Eigen]values Specified

            dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'Q',
                A(K + 1, 1), LDA, WORK(N + 1), IINFO);
          } else if (ITYPE == 5) {
            // Symmetric, eigenvalues specified

            dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, K, K, 'Q',
                A, LDA, WORK(N + 1), IINFO);
          } else if (ITYPE == 7) {
            // Diagonal, random eigenvalues

            dlatmr(
                N,
                N,
                'S',
                ISEED,
                'S',
                WORK,
                6,
                ONE,
                ONE,
                'T',
                'N',
                WORK(N + 1),
                1,
                ONE,
                WORK(2 * N + 1),
                1,
                ONE,
                'N',
                IDUMMA,
                0,
                0,
                ZERO,
                ANORM,
                'Q',
                A(K + 1, 1),
                LDA,
                IDUMMA,
                IINFO);
          } else if (ITYPE == 8) {
            // Symmetric, random eigenvalues

            dlatmr(
                N,
                N,
                'S',
                ISEED,
                'S',
                WORK,
                6,
                ONE,
                ONE,
                'T',
                'N',
                WORK(N + 1),
                1,
                ONE,
                WORK(2 * N + 1),
                1,
                ONE,
                'N',
                IDUMMA,
                K,
                K,
                ZERO,
                ANORM,
                'Q',
                A,
                LDA,
                IDUMMA,
                IINFO);
          } else if (ITYPE == 9) {
            // Positive definite, eigenvalues specified.

            dlatms(N, N, 'S', ISEED, 'P', WORK, IMODE, COND, ANORM, K, K, 'Q',
                A, LDA, WORK(N + 1), IINFO);
          } else if (ITYPE == 10) {
            // Positive definite tridiagonal, eigenvalues specified.

            if (N > 1) K = max(1, K);
            dlatms(N, N, 'S', ISEED, 'P', WORK, IMODE, COND, ANORM, 1, 1, 'Q',
                A(K, 1), LDA, WORK(N + 1), IINFO);
            for (I = 2; I <= N; I++) {
              TEMP1 =
                  A[K][I].abs() / sqrt((A[K + 1][I - 1] * A[K + 1][I]).abs());
              if (TEMP1 > HALF) {
                A[K][I] = HALF * sqrt((A[K + 1][I - 1] * A[K + 1][I]).abs());
              }
            }
          } else {
            IINFO.value = 1;
          }

          if (IINFO.value != 0) {
            print9999(NOUNIT, 'Generator', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            return;
          }
        }

        while (true) {
          // Call DSBTRD to compute S and U from upper triangle.

          dlacpy(' ', K + 1, N, A, LDA, WORK.asMatrix(LDA), LDA);

          NTEST = 1;
          dsbtrd('V', 'U', N, K, WORK.asMatrix(LDA), LDA, SD, SE, U, LDU,
              WORK(LDA * N + 1), IINFO);

          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSBTRD(U)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[1] = ULPINV;
              break;
            }
          }

          // Do tests 1 and 2

          dsbt21('Upper', N, K, 1, A, LDA, SD, SE, U, LDU, WORK, RESULT(1));

          // Convert A from Upper-Triangle-Only storage to
          // Lower-Triangle-Only storage.

          for (JC = 1; JC <= N; JC++) {
            for (JR = 0; JR <= min(K, N - JC); JR++) {
              A[JR + 1][JC] = A[K + 1 - JR][JC + JR];
            }
          }
          for (JC = N + 1 - K; JC <= N; JC++) {
            for (JR = min(K, N - JC) + 1; JR <= K; JR++) {
              A[JR + 1][JC] = ZERO;
            }
          }

          // Call DSBTRD to compute S and U from lower triangle

          dlacpy(' ', K + 1, N, A, LDA, WORK.asMatrix(LDA), LDA);

          NTEST = 3;
          dsbtrd('V', 'L', N, K, WORK.asMatrix(LDA), LDA, SD, SE, U, LDU,
              WORK(LDA * N + 1), IINFO);

          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSBTRD(L)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[3] = ULPINV;
              break;
            }
          }
          NTEST = 4;

          // Do tests 3 and 4

          dsbt21('Lower', N, K, 1, A, LDA, SD, SE, U, LDU, WORK, RESULT(3));

          // End of Loop -- Check for RESULT[j] > THRESH
        }
        NTESTT += NTEST;

        // Print out tests which fail.

        for (JR = 1; JR <= NTEST; JR++) {
          if (RESULT[JR] >= THRESH) {
            // If this is the first test to fail,
            // print a header to the data file.

            if (NERRS == 0) {
              NOUNIT.println(
                  '\n DSB -- Real Symmetric Banded Tridiagonal Reduction Routines');
              NOUNIT.println(' Matrix types (see DCHKSB for details): ');
              NOUNIT.println(
                  '\n Special Matrices:\n  1=Zero matrix.                          5=Diagonal: clustered entries.\n  2=Identity matrix.                      6=Diagonal: large, evenly spaced.\n  3=Diagonal: evenly spaced entries.      7=Diagonal: small, evenly spaced.\n  4=Diagonal: geometr. spaced entries.');
              NOUNIT.println(
                  ' Dense Symmetric Banded Matrices:\n  8=Evenly spaced eigenvals.             12=Small, evenly spaced eigenvals.\n  9=Geometrically spaced eigenvals.      13=Matrix with random O(1) entries.\n 10=Clustered eigenvalues.               14=Matrix with large random entries.\n 11=Large, evenly spaced eigenvals.      15=Matrix with small random entries.');
              NOUNIT.println(
                  '\n Tests performed:   (S is Tridiag,  U is orthogonal,\n${' ' * 20}\' means transpose.\n UPLO='
                  'U'
                  ':\n  1= | A - U S U\' | / ( |A| n ulp )       2= | I - U U\' | / ( n ulp )\n UPLO='
                  'L'
                  ':\n  3= | A - U S U\' | / ( |A| n ulp )       4= | I - U U\' | / ( n ulp )');
            }
            NERRS++;
            NOUNIT.println(
                ' N=${N.i5}, K=${K.i4}, seed=${IOLDSD.i4(4, ',')} type ${JTYPE.i2}, test(${JR.i2})=${RESULT[JR].g10_3}');
          }
        }
      }
    }
  }

  // Summary
  dlasum('DSB', NOUNIT, NERRS, NTESTT);
}

void print9999(
  final Nout NOUNIT,
  final String s,
  final int info,
  final int n,
  final int jtype,
  final Array<int> iseed,
) {
  NOUNIT.println(
      ' DCHKSB: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}
