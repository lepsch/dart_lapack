// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zbdsqr.dart';
import 'package:dart_lapack/src/zgebrd.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlaset.dart';
import 'package:dart_lapack/src/zungbr.dart';

import '../matgen/dlarnd.dart';
import '../matgen/zlatmr.dart';
import '../matgen/zlatms.dart';
import 'alasum.dart';
import 'common.dart';
import 'dlahd2.dart';
import 'dsvdch.dart';
import 'zbdt01.dart';
import 'zbdt02.dart';
import 'zbdt03.dart';
import 'zunt01.dart';

void zchkbd(
  final int NSIZES,
  final Array<int> MVAL_,
  final Array<int> NVAL_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final int NRHS,
  final Array<int> ISEED_,
  final double THRESH,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> BD_,
  final Array<double> BE_,
  final Array<double> S1_,
  final Array<double> S2_,
  final Matrix<Complex> X_,
  final int LDX,
  final Matrix<Complex> Y_,
  final Matrix<Complex> Z_,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Matrix<Complex> PT_,
  final int LDPT,
  final Matrix<Complex> U_,
  final Matrix<Complex> VT_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Nout NOUT,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final MVAL = MVAL_.having();
  final NVAL = NVAL_.having();
  final DOTYPE = DOTYPE_.having();
  final ISEED = ISEED_.having(length: 4);
  final A = A_.having(ld: LDA);
  final X = X_.having(ld: LDX);
  final Y = Y_.having(ld: LDX);
  final Z = Z_.having(ld: LDX);
  final Q = Q_.having(ld: LDQ);
  final PT = PT_.having(ld: LDPT);
  final U = U_.having(ld: LDPT);
  final VT = VT_.having(ld: LDPT);
  final WORK = WORK_.having();
  final BD = BD_.having();
  final BE = BE_.having();
  final S1 = S1_.having();
  final S2 = S2_.having();
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, HALF = 0.5;
  const MAXTYP = 16;
  bool BADMM, BADNN, BIDIAG = false;
  String UPLO;
  String PATH;
  int I,
      IMODE,
      ITYPE,
      J,
      JCOL,
      JSIZE,
      JTYPE,
      LOG2UI,
      M,
      MINWRK,
      MMAX,
      MNMAX,
      MNMIN = 0,
      MQ = 0,
      MTYPES,
      N,
      NFAIL,
      NMAX,
      NTEST;
  double AMNINV,
      ANORM = 0,
      COND,
      OVFL,
      RTOVFL,
      RTUNFL,
      TEMP1 = 0,
      TEMP2,
      ULP,
      ULPINV,
      UNFL;
  final IOLDSD = Array<int>(4), IWORK = Array<int>(1);
  final DUMMA = Array<double>(1), RESULT = Array<double>(14);
  final IINFO = Box(0);
  final KTYPE = Array.fromList([
    1, 2, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 9, 9, 9, 10 //
  ]);
  final KMAGN = Array.fromList([
    1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 0 //
  ]);
  final KMODE = Array.fromList([
    0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 0 //
  ]);

  // Check for errors

  INFO.value = 0;

  BADMM = false;
  BADNN = false;
  MMAX = 1;
  NMAX = 1;
  MNMAX = 1;
  MINWRK = 1;
  for (J = 1; J <= NSIZES; J++) {
    MMAX = max(MMAX, MVAL[J]);
    if (MVAL[J] < 0) BADMM = true;
    NMAX = max(NMAX, NVAL[J]);
    if (NVAL[J] < 0) BADNN = true;
    MNMAX = max(MNMAX, min(MVAL[J], NVAL[J]));
    MINWRK = max(
        MINWRK,
        max(
            3 * (MVAL[J] + NVAL[J]),
            MVAL[J] * (MVAL[J] + max(MVAL[J], max(NVAL[J], NRHS)) + 1).toInt() +
                NVAL[J] * min(NVAL[J], MVAL[J])));
  }

  // Check for errors

  if (NSIZES < 0) {
    INFO.value = -1;
  } else if (BADMM) {
    INFO.value = -2;
  } else if (BADNN) {
    INFO.value = -3;
  } else if (NTYPES < 0) {
    INFO.value = -4;
  } else if (NRHS < 0) {
    INFO.value = -6;
  } else if (LDA < MMAX) {
    INFO.value = -11;
  } else if (LDX < MMAX) {
    INFO.value = -17;
  } else if (LDQ < MMAX) {
    INFO.value = -21;
  } else if (LDPT < MNMAX) {
    INFO.value = -23;
  } else if (MINWRK > LWORK) {
    INFO.value = -27;
  }

  if (INFO.value != 0) {
    xerbla('ZCHKBD', -INFO.value);
    return;
  }

  // Initialize constants

  PATH = '${'Zomplex precision'[0]}BD';
  NFAIL = 0;
  NTEST = 0;
  UNFL = dlamch('Safe minimum');
  OVFL = dlamch('Overflow');
  ULP = dlamch('Precision');
  ULPINV = ONE / ULP;
  LOG2UI = log(ULPINV) ~/ log(TWO);
  RTUNFL = sqrt(UNFL);
  RTOVFL = sqrt(OVFL);
  infoc.INFOT = 0;

  // Loop over sizes, types

  for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) {
    M = MVAL[JSIZE];
    N = NVAL[JSIZE];
    MNMIN = min(M, N);
    AMNINV = ONE / max(M, max(N, 1));

    if (NSIZES != 1) {
      MTYPES = min(MAXTYP, NTYPES);
    } else {
      MTYPES = min(MAXTYP + 1, NTYPES);
    }

    for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
      if (!DOTYPE[JTYPE]) continue;

      for (J = 1; J <= 4; J++) {
        IOLDSD[J] = ISEED[J];
      }

      for (J = 1; J <= 14; J++) {
        RESULT[J] = -ONE;
      }

      UPLO = ' ';

      // Compute "A"

      // Control parameters:

      // KMAGN  KMODE        KTYPE
      //    =1  O(1)   clustered 1  zero
      //    =2  large  clustered 2  identity
      //    =3  small  exponential  (none)
      //    =4         arithmetic   diagonal, (w/ eigenvalues)
      //    =5         random       symmetric, w/ eigenvalues
      //    =6                      nonsymmetric, w/ singular values
      //    =7                      random diagonal
      //    =8                      random symmetric
      //    =9                      random nonsymmetric
      //    =10                     random bidiagonal (log. distrib.)

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
        IINFO.value = 0;
        COND = ULPINV;

        BIDIAG = false;
        if (ITYPE == 1) {
          // Zero matrix

          IINFO.value = 0;
        } else if (ITYPE == 2) {
          // Identity

          for (JCOL = 1; JCOL <= MNMIN; JCOL++) {
            A[JCOL][JCOL] = ANORM.toComplex();
          }
        } else if (ITYPE == 4) {
          // Diagonal Matrix, [Eigen]values Specified

          zlatms(MNMIN, MNMIN, 'S', ISEED, 'N', RWORK, IMODE, COND, ANORM, 0, 0,
              'N', A, LDA, WORK, IINFO);
        } else if (ITYPE == 5) {
          // Symmetric, eigenvalues specified

          zlatms(MNMIN, MNMIN, 'S', ISEED, 'S', RWORK, IMODE, COND, ANORM, M, N,
              'N', A, LDA, WORK, IINFO);
        } else if (ITYPE == 6) {
          // Nonsymmetric, singular values specified

          zlatms(M, N, 'S', ISEED, 'N', RWORK, IMODE, COND, ANORM, M, N, 'N', A,
              LDA, WORK, IINFO);
        } else if (ITYPE == 7) {
          // Diagonal, random entries

          zlatmr(
              MNMIN,
              MNMIN,
              'S',
              ISEED,
              'N',
              WORK,
              6,
              ONE,
              Complex.one,
              'T',
              'N',
              WORK(MNMIN + 1),
              1,
              ONE,
              WORK(2 * MNMIN + 1),
              1,
              ONE,
              'N',
              IWORK,
              0,
              0,
              ZERO,
              ANORM,
              'NO',
              A,
              LDA,
              IWORK,
              IINFO);
        } else if (ITYPE == 8) {
          // Symmetric, random entries

          zlatmr(
              MNMIN,
              MNMIN,
              'S',
              ISEED,
              'S',
              WORK,
              6,
              ONE,
              Complex.one,
              'T',
              'N',
              WORK(MNMIN + 1),
              1,
              ONE,
              WORK(M + MNMIN + 1),
              1,
              ONE,
              'N',
              IWORK,
              M,
              N,
              ZERO,
              ANORM,
              'NO',
              A,
              LDA,
              IWORK,
              IINFO);
        } else if (ITYPE == 9) {
          // Nonsymmetric, random entries

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
              WORK(MNMIN + 1),
              1,
              ONE,
              WORK(M + MNMIN + 1),
              1,
              ONE,
              'N',
              IWORK,
              M,
              N,
              ZERO,
              ANORM,
              'NO',
              A,
              LDA,
              IWORK,
              IINFO);
        } else if (ITYPE == 10) {
          // Bidiagonal, random entries

          TEMP1 = -TWO * log(ULP);
          for (J = 1; J <= MNMIN; J++) {
            BD[J] = exp(TEMP1 * dlarnd(2, ISEED));
            if (J < MNMIN) BE[J] = exp(TEMP1 * dlarnd(2, ISEED));
          }

          IINFO.value = 0;
          BIDIAG = true;
          if (M >= N) {
            UPLO = 'U';
          } else {
            UPLO = 'L';
          }
        } else {
          IINFO.value = 1;
        }

        if (IINFO.value == 0) {
          // Generate Right-Hand Side

          if (BIDIAG) {
            zlatmr(
                MNMIN,
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
                WORK(MNMIN + 1),
                1,
                ONE,
                WORK(2 * MNMIN + 1),
                1,
                ONE,
                'N',
                IWORK,
                MNMIN,
                NRHS,
                ZERO,
                ONE,
                'NO',
                Y,
                LDX,
                IWORK,
                IINFO);
          } else {
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
                IWORK,
                M,
                NRHS,
                ZERO,
                ONE,
                'NO',
                X,
                LDX,
                IWORK,
                IINFO);
          }
        }

        // Error Exit

        if (IINFO.value != 0) {
          _print9998(NOUT, 'Generator', IINFO.value, M, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }
      }

      // Call ZGEBRD and ZUNGBR to compute B, Q, and P, do tests.

      if (!BIDIAG) {
        // Compute transformations to reduce A to bidiagonal form:
        // B := Q' * A * P.

        zlacpy(' ', M, N, A, LDA, Q, LDQ);
        zgebrd(M, N, Q, LDQ, BD, BE, WORK, WORK(MNMIN + 1), WORK(2 * MNMIN + 1),
            LWORK - 2 * MNMIN, IINFO);

        // Check error code from ZGEBRD.

        if (IINFO.value != 0) {
          _print9998(NOUT, 'ZGEBRD', IINFO.value, M, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }

        zlacpy(' ', M, N, Q, LDQ, PT, LDPT);
        if (M >= N) {
          UPLO = 'U';
        } else {
          UPLO = 'L';
        }

        // Generate Q

        MQ = M;
        if (NRHS <= 0) MQ = MNMIN;
        zungbr('Q', M, MQ, N, Q, LDQ, WORK, WORK(2 * MNMIN + 1),
            LWORK - 2 * MNMIN, IINFO);

        // Check error code from ZUNGBR.

        if (IINFO.value != 0) {
          _print9998(NOUT, 'ZUNGBR(Q)', IINFO.value, M, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }

        // Generate P'

        zungbr('P', MNMIN, N, M, PT, LDPT, WORK(MNMIN + 1), WORK(2 * MNMIN + 1),
            LWORK - 2 * MNMIN, IINFO);

        // Check error code from ZUNGBR.

        if (IINFO.value != 0) {
          _print9998(NOUT, 'ZUNGBR(P)', IINFO.value, M, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }

        // Apply Q' to an M by NRHS matrix X:  Y := Q' * X.

        zgemm('Conjugate transpose', 'No transpose', M, NRHS, M, Complex.one, Q,
            LDQ, X, LDX, Complex.zero, Y, LDX);

        // Test 1:  Check the decomposition A := Q * B * PT
        //      2:  Check the orthogonality of Q
        //      3:  Check the orthogonality of PT

        zbdt01(M, N, 1, A, LDA, Q, LDQ, BD, BE, PT, LDPT, WORK, RWORK,
            RESULT.box(1));
        zunt01('Columns', M, MQ, Q, LDQ, WORK, LWORK, RWORK, RESULT.box(2));
        zunt01('Rows', MNMIN, N, PT, LDPT, WORK, LWORK, RWORK, RESULT.box(3));
      }

      // Use ZBDSQR to form the SVD of the bidiagonal matrix B:
      // B := U * S1 * VT, and compute Z = U' * Y.

      dcopy(MNMIN, BD, 1, S1, 1);
      if (MNMIN > 0) dcopy(MNMIN - 1, BE, 1, RWORK, 1);
      zlacpy(' ', M, NRHS, Y, LDX, Z, LDX);
      zlaset('Full', MNMIN, MNMIN, Complex.zero, Complex.one, U, LDPT);
      zlaset('Full', MNMIN, MNMIN, Complex.zero, Complex.one, VT, LDPT);

      zbdsqr(UPLO, MNMIN, MNMIN, MNMIN, NRHS, S1, RWORK, VT, LDPT, U, LDPT, Z,
          LDX, RWORK(MNMIN + 1), IINFO);

      // Check error code from ZBDSQR.

      if (IINFO.value != 0) {
        _print9998(NOUT, 'ZBDSQR(vects)', IINFO.value, M, N, JTYPE, IOLDSD);
        INFO.value = (IINFO.value).abs();
        if (IINFO.value < 0) {
          return;
        }

        RESULT[4] = ULPINV;
      } else {
        // Use ZBDSQR to compute only the singular values of the
        // bidiagonal matrix B;  U, VT, and Z should not be modified.

        dcopy(MNMIN, BD, 1, S2, 1);
        if (MNMIN > 0) dcopy(MNMIN - 1, BE, 1, RWORK, 1);

        zbdsqr(UPLO, MNMIN, 0, 0, 0, S2, RWORK, VT, LDPT, U, LDPT, Z, LDX,
            RWORK(MNMIN + 1), IINFO);

        // Check error code from ZBDSQR.

        if (IINFO.value != 0) {
          _print9998(NOUT, 'ZBDSQR(values)', IINFO.value, M, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          }

          RESULT[9] = ULPINV;
        } else {
          // Test 4:  Check the decomposition B := U * S1 * VT
          //      5:  Check the computation Z := U' * Y
          //      6:  Check the orthogonality of U
          //      7:  Check the orthogonality of VT

          zbdt03(UPLO, MNMIN, 1, BD, BE, U, LDPT, S1, VT, LDPT, WORK,
              RESULT.box(4));
          zbdt02(
              MNMIN, NRHS, Y, LDX, Z, LDX, U, LDPT, WORK, RWORK, RESULT.box(5));
          zunt01('Columns', MNMIN, MNMIN, U, LDPT, WORK, LWORK, RWORK,
              RESULT.box(6));
          zunt01('Rows', MNMIN, MNMIN, VT, LDPT, WORK, LWORK, RWORK,
              RESULT.box(7));

          // Test 8:  Check that the singular values are sorted in
          //          non-increasing order and are non-negative

          RESULT[8] = ZERO;
          for (I = 1; I <= MNMIN - 1; I++) {
            if (S1[I] < S1[I + 1]) RESULT[8] = ULPINV;
            if (S1[I] < ZERO) RESULT[8] = ULPINV;
          }
          if (MNMIN >= 1) {
            if (S1[MNMIN] < ZERO) RESULT[8] = ULPINV;
          }

          // Test 9:  Compare ZBDSQR with and without singular vectors

          TEMP2 = ZERO;

          for (J = 1; J <= MNMIN; J++) {
            TEMP1 = (S1[J] - S2[J]).abs() /
                max(sqrt(UNFL) * max(S1[1], ONE),
                    ULP * max(S1[J].abs(), S2[J].abs()));
            TEMP2 = max(TEMP1, TEMP2);
          }

          RESULT[9] = TEMP2;

          // Test 10:  Sturm sequence test of singular values
          //           Go up by factors of two until it succeeds

          TEMP1 = THRESH * (HALF - ULP);

          for (J = 0; J <= LOG2UI; J++) {
            dsvdch(MNMIN, BD, BE, S1, TEMP1, IINFO);
            if (IINFO.value == 0) break;
            TEMP1 *= TWO;
          }
          RESULT[10] = TEMP1;

          // Use ZBDSQR to form the decomposition A := (QU) S (VT PT)
          // from the bidiagonal form A := Q B PT.

          if (!BIDIAG) {
            dcopy(MNMIN, BD, 1, S2, 1);
            if (MNMIN > 0) dcopy(MNMIN - 1, BE, 1, RWORK, 1);

            zbdsqr(UPLO, MNMIN, N, M, NRHS, S2, RWORK, PT, LDPT, Q, LDQ, Y, LDX,
                RWORK(MNMIN + 1), IINFO);

            // Test 11:  Check the decomposition A := Q*U * S2 * VT*PT
            //      12:  Check the computation Z := U' * Q' * X
            //      13:  Check the orthogonality of Q*U
            //      14:  Check the orthogonality of VT*PT

            zbdt01(M, N, 0, A, LDA, Q, LDQ, S2, DUMMA, PT, LDPT, WORK, RWORK,
                RESULT.box(11));
            zbdt02(
                M, NRHS, X, LDX, Y, LDX, Q, LDQ, WORK, RWORK, RESULT.box(12));
            zunt01(
                'Columns', M, MQ, Q, LDQ, WORK, LWORK, RWORK, RESULT.box(13));
            zunt01(
                'Rows', MNMIN, N, PT, LDPT, WORK, LWORK, RWORK, RESULT.box(14));
          }
        }
      }

      // End of Loop -- Check for RESULT(j) > THRESH

      for (J = 1; J <= 14; J++) {
        if (RESULT[J] >= THRESH) {
          if (NFAIL == 0) dlahd2(NOUT, PATH);
          NOUT.println(
              ' M=${M.i5}, N=${N.i5}, type ${JTYPE.i2}, seed=${IOLDSD.i4(4, ',')} test(${J.i2})=${RESULT[J].g11_4}');
          NFAIL++;
        }
      }
      if (!BIDIAG) {
        NTEST += 14;
      } else {
        NTEST += 5;
      }
    }
  }

  // Summary

  alasum(PATH, NOUT, NFAIL, NTEST, 0);
}

void _print9998(
    Nout nout, String s, int info, int m, int n, int jtype, Array<int> iseed) {
  nout.println(
      ' ZCHKBD: $s returned INFO=${info.i6}.\n${' ' * 9}M=${m.i6}, N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}
