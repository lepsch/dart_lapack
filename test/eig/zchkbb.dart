// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zgbbrd.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlaset.dart';

import '../matgen/zlatmr.dart';
import '../matgen/zlatms.dart';
import 'dlahd2.dart';
import 'dlasum.dart';
import 'zbdt01.dart';
import 'zbdt02.dart';
import 'zunt01.dart';

void zchkbb(
  final int NSIZES,
  final Array<int> MVAL_,
  final Array<int> NVAL_,
  final int NWDTHS,
  final Array<int> KK_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final int NRHS,
  final Array<int> ISEED_,
  final double THRESH,
  final Nout NOUNIT,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Array<double> BD_,
  final Array<double> BE_,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Matrix<Complex> P_,
  final int LDP,
  final Matrix<Complex> C_,
  final int LDC,
  final Matrix<Complex> CC_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<double> RESULT_,
  final Box<int> INFO,
) {
  final MVAL = MVAL_.having();
  final NVAL = NVAL_.having();
  final KK = KK_.having();
  final DOTYPE = DOTYPE_.having();
  final ISEED = ISEED_.having(length: 4);
  final A = A_.having(ld: LDA);
  final AB = AB_.having(ld: LDAB);
  final Q = Q_.having(ld: LDQ);
  final P = P_.having(ld: LDP);
  final C = C_.having(ld: LDC);
  final CC = CC_.having(ld: LDC);
  final BD = BD_.having();
  final BE = BE_.having();
  final RESULT = RESULT_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();

// -- LAPACK test routine (input) --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  const MAXTYP = 15;
  bool BADMM, BADNN, BADNNB;
  int I,
      IMODE,
      ITYPE,
      J,
      JCOL,
      JR,
      JSIZE,
      JTYPE,
      JWIDTH,
      K,
      KL,
      KMAX,
      KU,
      M,
      MMAX,
      MNMAX,
      // MNMIN,
      MTYPES,
      N,
      NERRS,
      NMAX,
      NTEST,
      NTESTT;
  double AMNINV, ANORM = 0, COND, OVFL, RTOVFL, RTUNFL, ULP, ULPINV, UNFL;
  final IINFO = Box(0);
  final IDUMMA = Array<int>(1), IOLDSD = Array<int>(4);
  final KTYPE = Array.fromList([
    1, 2, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 9, 9, 9, //
  ]);
  final KMAGN = Array.fromList([
    1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3 //
  ]);
  final KMODE = Array.fromList([
    0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0 //
  ]);

  // Check for errors

  NTESTT = 0;
  INFO.value = 0;

  // Important constants

  BADMM = false;
  BADNN = false;
  MMAX = 1;
  NMAX = 1;
  MNMAX = 1;
  for (J = 1; J <= NSIZES; J++) {
    MMAX = max(MMAX, MVAL[J]);
    if (MVAL[J] < 0) BADMM = true;
    NMAX = max(NMAX, NVAL[J]);
    if (NVAL[J] < 0) BADNN = true;
    MNMAX = max(MNMAX, min(MVAL[J], NVAL[J]));
  }

  BADNNB = false;
  KMAX = 0;
  for (J = 1; J <= NWDTHS; J++) {
    KMAX = max(KMAX, KK[J]);
    if (KK[J] < 0) BADNNB = true;
  }

  // Check for errors

  if (NSIZES < 0) {
    INFO.value = -1;
  } else if (BADMM) {
    INFO.value = -2;
  } else if (BADNN) {
    INFO.value = -3;
  } else if (NWDTHS < 0) {
    INFO.value = -4;
  } else if (BADNNB) {
    INFO.value = -5;
  } else if (NTYPES < 0) {
    INFO.value = -6;
  } else if (NRHS < 0) {
    INFO.value = -8;
  } else if (LDA < NMAX) {
    INFO.value = -13;
  } else if (LDAB < 2 * KMAX + 1) {
    INFO.value = -15;
  } else if (LDQ < NMAX) {
    INFO.value = -19;
  } else if (LDP < NMAX) {
    INFO.value = -21;
  } else if (LDC < NMAX) {
    INFO.value = -23;
  } else if ((max(LDA, NMAX) + 1) * NMAX > LWORK) {
    INFO.value = -26;
  }

  if (INFO.value != 0) {
    xerbla('ZCHKBB', -INFO.value);
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

  // Loop over sizes, widths, types

  NERRS = 0;

  for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) {
    M = MVAL[JSIZE];
    N = NVAL[JSIZE];
    //  MNMIN = min( M, N );
    AMNINV = ONE / max(1, max(M, N));

    for (JWIDTH = 1; JWIDTH <= NWDTHS; JWIDTH++) {
      K = KK[JWIDTH];
      if (K >= M && K >= N) continue;
      KL = max(0, min(M - 1, K));
      KU = max(0, min(N - 1, K));

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

        // Control parameters:
        //
        //     KMAGN  KMODE        KTYPE
        //        =1  O(1)   clustered 1  zero
        //        =2  large  clustered 2  identity
        //        =3  small  exponential  (none)
        //        =4         arithmetic   diagonal, (w/ singular values)
        //        =5         random log   (none)
        //        =6         random       nonhermitian, w/ singular values
        //        =7                      (none)
        //        =8                      (none)
        //        =9                      random nonhermitian

        if (MTYPES <= MAXTYP) {
          ITYPE = KTYPE[JTYPE];
          IMODE = KMODE[JTYPE];

          // Compute norm

          switch (KMAGN[JTYPE]) {
            case 1:
              ANORM = ONE;
              break;

            case 2:
              ANORM = (RTOVFL * ULP) * AMNINV;
              break;

            case 3:
              ANORM = RTUNFL * max(M, N) * ULPINV;
              break;
          }

          zlaset('Full', LDA, N, Complex.zero, Complex.zero, A, LDA);
          zlaset('Full', LDAB, N, Complex.zero, Complex.zero, AB, LDAB);
          IINFO.value = 0;
          COND = ULPINV;

          // Special Matrices -- Identity & Jordan block

          //    Zero

          if (ITYPE == 1) {
            IINFO.value = 0;
          } else if (ITYPE == 2) {
            // Identity

            for (JCOL = 1; JCOL <= N; JCOL++) {
              A[JCOL][JCOL] = ANORM.toComplex();
            }
          } else if (ITYPE == 4) {
            // Diagonal Matrix, singular values specified

            zlatms(M, N, 'S', ISEED, 'N', RWORK, IMODE, COND, ANORM, 0, 0, 'N',
                A, LDA, WORK, IINFO);
          } else if (ITYPE == 6) {
            // Nonhermitian, singular values specified

            zlatms(M, N, 'S', ISEED, 'N', RWORK, IMODE, COND, ANORM, KL, KU,
                'N', A, LDA, WORK, IINFO);
          } else if (ITYPE == 9) {
            // Nonhermitian, random entries

            zlatmr(
                M,
                N,
                'S',
                ISEED,
                'N',
                WORK,
                6,
                ONE,
                Complex.one,
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
                KL,
                KU,
                ZERO,
                ANORM,
                'N',
                A,
                LDA,
                IDUMMA,
                IINFO);
          } else {
            IINFO.value = 1;
          }

          // Generate Right-Hand Side

          zlatmr(
              M,
              NRHS,
              'S',
              ISEED,
              'N',
              WORK,
              6,
              ONE,
              Complex.one,
              'T',
              'N',
              WORK(M + 1),
              1,
              ONE,
              WORK(2 * M + 1),
              1,
              ONE,
              'N',
              IDUMMA,
              M,
              NRHS,
              ZERO,
              ONE,
              'NO',
              C,
              LDC,
              IDUMMA,
              IINFO);

          if (IINFO.value != 0) {
            _print9999(
                NOUNIT, 'Generator', IINFO.value, M, N, K, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            return;
          }
        }

        // Copy A to band storage.

        for (J = 1; J <= N; J++) {
          for (I = max(1, J - KU); I <= min(M, J + KL); I++) {
            AB[KU + 1 + I - J][J] = A[I][J];
          }
        }

        // Copy C

        zlacpy('Full', M, NRHS, C, LDC, CC, LDC);

        // Call ZGBBRD to compute B, Q and P, and to update C.

        zgbbrd('B', M, N, NRHS, KL, KU, AB, LDAB, BD, BE, Q, LDQ, P, LDP, CC,
            LDC, WORK, RWORK, IINFO);

        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZGBBRD', IINFO.value, M, N, K, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          }
          RESULT[1] = ULPINV;
        } else {
          // Test 1:  Check the decomposition A := Q * B * P'
          //      2:  Check the orthogonality of Q
          //      3:  Check the orthogonality of P
          //      4:  Check the computation of Q' * C

          zbdt01(M, N, -1, A, LDA, Q, LDQ, BD, BE, P, LDP, WORK, RWORK,
              RESULT.box(1));
          zunt01('Columns', M, M, Q, LDQ, WORK, LWORK, RWORK, RESULT.box(2));
          zunt01('Rows', N, N, P, LDP, WORK, LWORK, RWORK, RESULT.box(3));
          zbdt02(M, NRHS, C, LDC, CC, LDC, Q, LDQ, WORK, RWORK, RESULT.box(4));

          // End of Loop -- Check for RESULT(j) > THRESH

          NTEST = 4;
        }
        NTESTT += NTEST;

        // Print out tests which fail.

        for (JR = 1; JR <= NTEST; JR++) {
          if (RESULT[JR] >= THRESH) {
            if (NERRS == 0) dlahd2(NOUNIT, 'ZBB');
            NERRS++;
            NOUNIT.println(
                ' M =${M.i4} N=${N.i4}, K=${K.i3}, seed=${IOLDSD.i4(4, ',')} type ${JTYPE.i2}, test(${JR.i2})=${RESULT[JR].g10_3}');
          }
        }
      }
    }

    // Summary

    dlasum('ZBB', NOUNIT, NERRS, NTESTT);
  }
}

void _print9999(Nout nout, String s, int info, int m, int n, int k, int jtype,
    Array<int> iseed) {
  nout.println(
      ' ZCHKBB: $s returned INFO=${info.i5}.\n${' ' * 9}M=${m.i5} N=${n.i5} K=${k.i5}, JTYPE=${jtype.i5}, ISEED=(${iseed.i5(4, ',')})');
}
