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
import 'package:dart_lapack/src/zlaset.dart';

import '../matgen/zlatme.dart';
import '../matgen/zlatmr.dart';
import '../matgen/zlatms.dart';
import 'dlasum.dart';
import 'zget24.dart';

Future<void> zdrvsx(
  final int NSIZES,
  final Array<int> NN_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final Array<int> ISEED_,
  final double THRESH,
  final Nin NIUNIT,
  final Nout NOUNIT,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> H_,
  final Matrix<Complex> HT_,
  final Array<Complex> W_,
  final Array<Complex> WT_,
  final Array<Complex> WTMP_,
  final Matrix<Complex> VS_,
  final int LDVS,
  final Matrix<Complex> VS1_,
  final Array<double> RESULT_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<bool> BWORK_,
  final Box<int> INFO,
) async {
  final NN = NN_.having();
  final DOTYPE = DOTYPE_.having();
  final ISEED = ISEED_.having(length: 4);
  final A = A_.having(ld: LDA);
  final H = H_.having(ld: LDA);
  final HT = HT_.having(ld: LDA);
  final VS = VS_.having(ld: LDVS);
  final VS1 = VS1_.having(ld: LDVS);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final BWORK = BWORK_.having();
  final W = W_.having();
  final WT = WT_.having();
  final WTMP = WTMP_.having();
  final RESULT = RESULT_.having(length: 17);
  const ZERO = 0.0, ONE = 1.0;
  const MAXTYP = 21;
  bool BADNN;
  String PATH;
  int IMODE,
      ISRT,
      ITYPE,
      IWK,
      J,
      JCOL,
      JSIZE,
      JTYPE,
      MTYPES,
      N,
      NERRS = 0,
      NFAIL,
      NMAX,
      NNWORK,
      NSLCT = 0,
      NTEST,
      NTESTF,
      NTESTT;
  double ANORM = 0,
      COND,
      CONDS,
      OVFL,
      RCDEIN = 0,
      RCDVIN = 0,
      RTULP,
      RTULPI,
      ULP,
      ULPINV,
      UNFL;
  final IDUMMA = Array<int>(1), IOLDSD = Array<int>(4), ISLCT = Array<int>(20);
  final KTYPE = Array.fromList(
      [1, 2, 3, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9]);
  final KMAGN = Array.fromList(
      [1, 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 1, 2, 3]);
  final KMODE = Array.fromList(
      [0, 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1]);
  final KCONDS = Array.fromList(
      [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0, 0, 0]);
  final IINFO = Box(0);

  PATH = '${'Zomplex precision'[0]}SX';

  // Check for errors

  NTESTT = 0;
  NTESTF = 0;
  INFO.value = 0;

  // Important constants

  BADNN = false;

  // 8 is the largest dimension in the input file of precomputed
  // problems

  NMAX = 8;
  for (J = 1; J <= NSIZES; J++) {
    NMAX = max(NMAX, NN[J]);
    if (NN[J] < 0) BADNN = true;
  }

  // Check for errors

  if (NSIZES < 0) {
    INFO.value = -1;
  } else if (BADNN) {
    INFO.value = -2;
  } else if (NTYPES < 0) {
    INFO.value = -3;
  } else if (THRESH < ZERO) {
    INFO.value = -6;
    // } else if ( NIUNIT <= 0 ) {
    //    INFO = -7;
    // } else if ( NOUNIT <= 0 ) {
    //    INFO = -8;
  } else if (LDA < 1 || LDA < NMAX) {
    INFO.value = -10;
  } else if (LDVS < 1 || LDVS < NMAX) {
    INFO.value = -20;
  } else if (max(3 * NMAX, 2 * pow(NMAX, 2)) > LWORK) {
    INFO.value = -24;
  }

  if (INFO.value != 0) {
    xerbla('ZDRVSX', -INFO.value);
    return;
  }

  // If nothing to do check on NIUNIT

  if (NSIZES != 0 && NTYPES != 0) {
    // More Important constants

    UNFL = dlamch('Safe minimum');
    OVFL = ONE / UNFL;
    ULP = dlamch('Precision');
    ULPINV = ONE / ULP;
    RTULP = sqrt(ULP);
    RTULPI = ONE / RTULP;

    // Loop over sizes, types

    NERRS = 0;

    for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) {
      N = NN[JSIZE];
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

        // Compute "A"

        // Control parameters:

        //     KMAGN  KCONDS  KMODE        KTYPE
        // =1  O(1)   1       clustered 1  zero
        // =2  large  large   clustered 2  identity
        // =3  small          exponential  Jordan
        // =4                 arithmetic   diagonal, (w/ eigenvalues)
        // =5                 random log   symmetric, w/ eigenvalues
        // =6                 random       general, w/ eigenvalues
        // =7                              random diagonal
        // =8                              random symmetric
        // =9                              random general
        // =10                             random triangular

        if (MTYPES <= MAXTYP) {
          ITYPE = KTYPE[JTYPE];
          IMODE = KMODE[JTYPE];

          // Compute norm

          switch (KMAGN[JTYPE]) {
            case 1:
              ANORM = ONE;
              break;

            case 2:
              ANORM = OVFL * ULP;
              break;

            case 3:
              ANORM = UNFL * ULPINV;
              break;
          }

          zlaset('Full', LDA, N, Complex.zero, Complex.zero, A, LDA);
          IINFO.value = 0;
          COND = ULPINV;

          // Special Matrices -- Identity & Jordan block

          if (ITYPE == 1) {
            // Zero

            IINFO.value = 0;
          } else if (ITYPE == 2) {
            // Identity

            for (JCOL = 1; JCOL <= N; JCOL++) {
              A[JCOL][JCOL] = ANORM.toComplex();
            }
          } else if (ITYPE == 3) {
            // Jordan Block

            for (JCOL = 1; JCOL <= N; JCOL++) {
              A[JCOL][JCOL] = ANORM.toComplex();
              if (JCOL > 1) A[JCOL][JCOL - 1] = Complex.one;
            }
          } else if (ITYPE == 4) {
            // Diagonal Matrix, [Eigen]values Specified

            zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, 0, 0, 'N',
                A, LDA, WORK(N + 1), IINFO);
          } else if (ITYPE == 5) {
            // Symmetric, eigenvalues specified

            zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, N, N, 'N',
                A, LDA, WORK(N + 1), IINFO);
          } else if (ITYPE == 6) {
            // General, eigenvalues specified

            if (KCONDS[JTYPE] == 1) {
              CONDS = ONE;
            } else if (KCONDS[JTYPE] == 2) {
              CONDS = RTULPI;
            } else {
              CONDS = ZERO;
            }

            zlatme(N, 'D', ISEED, WORK, IMODE, COND, Complex.one, 'T', 'T', 'T',
                RWORK, 4, CONDS, N, N, ANORM, A, LDA, WORK(2 * N + 1), IINFO);
          } else if (ITYPE == 7) {
            // Diagonal, random eigenvalues

            zlatmr(
                N,
                N,
                'D',
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
                0,
                0,
                ZERO,
                ANORM,
                'NO',
                A,
                LDA,
                IDUMMA,
                IINFO);
          } else if (ITYPE == 8) {
            // Symmetric, random eigenvalues

            zlatmr(
                N,
                N,
                'D',
                ISEED,
                'H',
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
                N,
                N,
                ZERO,
                ANORM,
                'NO',
                A,
                LDA,
                IDUMMA,
                IINFO);
          } else if (ITYPE == 9) {
            // General, random eigenvalues

            zlatmr(
                N,
                N,
                'D',
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
                N,
                N,
                ZERO,
                ANORM,
                'NO',
                A,
                LDA,
                IDUMMA,
                IINFO);
            if (N >= 4) {
              zlaset('Full', 2, N, Complex.zero, Complex.zero, A, LDA);
              zlaset(
                  'Full', N - 3, 1, Complex.zero, Complex.zero, A(3, 1), LDA);
              zlaset('Full', N - 3, 2, Complex.zero, Complex.zero, A(3, N - 1),
                  LDA);
              zlaset('Full', 1, N, Complex.zero, Complex.zero, A(N, 1), LDA);
            }
          } else if (ITYPE == 10) {
            // Triangular, random eigenvalues

            zlatmr(
                N,
                N,
                'D',
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
                N,
                0,
                ZERO,
                ANORM,
                'NO',
                A,
                LDA,
                IDUMMA,
                IINFO);
          } else {
            IINFO.value = 1;
          }

          if (IINFO.value != 0) {
            NOUNIT.println(
                ' ZDRVSX: Generator returned INFO.value=${IINFO.value.i6}.\n${' ' * 9}N=${N.i6}, JTYPE=${JTYPE.i6}, ISEED=(${IOLDSD.i5(4, ',')})');
            INFO.value = (IINFO.value).abs();
            return;
          }
        }

        // Test for minimal and generous workspace

        for (IWK = 1; IWK <= 2; IWK++) {
          if (IWK == 1) {
            NNWORK = 2 * N;
          } else {
            NNWORK = max(2 * N, N * (N + 1) ~/ 2);
          }
          NNWORK = max(NNWORK, 1);

          zget24(
              false,
              JTYPE,
              THRESH,
              IOLDSD,
              NOUNIT,
              N,
              A,
              LDA,
              H,
              HT,
              W,
              WT,
              WTMP,
              VS,
              LDVS,
              VS1,
              RCDEIN,
              RCDVIN,
              NSLCT,
              ISLCT,
              0,
              RESULT,
              WORK,
              NNWORK,
              RWORK,
              BWORK,
              INFO);

          // Check for RESULT(j) > THRESH

          NTEST = 0;
          NFAIL = 0;
          for (J = 1; J <= 15; J++) {
            if (RESULT[J] >= ZERO) NTEST++;
            if (RESULT[J] >= THRESH) NFAIL++;
          }

          if (NFAIL > 0) NTESTF++;
          if (NTESTF == 1) {
            _printFirst(NOUNIT, PATH, THRESH);
            NTESTF = 2;
          }

          for (J = 1; J <= 15; J++) {
            if (RESULT[J] >= THRESH) {
              NOUNIT.println(
                  ' N=${N.i5}, IWK=${IWK.i2}, seed=${IOLDSD.i4(4, ',')} type ${JTYPE.i2}, test(${J.i2})=${RESULT[J].g10_3}');
            }
          }

          NERRS += NFAIL;
          NTESTT += NTEST;
        }
      }
    }
  }

  // Read in data from file to check accuracy of condition estimation
  // Read input data until N=0

  JTYPE = 0;
  try {
    while (true) {
      (N, NSLCT, ISRT) = await NIUNIT.readInt3();
      if (N == 0) break;

      JTYPE++;
      ISEED[1] = JTYPE;
      await NIUNIT.readArray(ISLCT, NSLCT);
      await NIUNIT.readMatrix(A, N, N);
      (RCDEIN, RCDVIN) = await NIUNIT.readDouble2();

      zget24(
          true,
          22,
          THRESH,
          ISEED,
          NOUNIT,
          N,
          A,
          LDA,
          H,
          HT,
          W,
          WT,
          WTMP,
          VS,
          LDVS,
          VS1,
          RCDEIN,
          RCDVIN,
          NSLCT,
          ISLCT,
          ISRT,
          RESULT,
          WORK,
          LWORK,
          RWORK,
          BWORK,
          INFO);

      // Check for RESULT(j) > THRESH

      NTEST = 0;
      var NFAIL = 0;
      for (J = 1; J <= 17; J++) {
        if (RESULT[J] >= ZERO) NTEST++;
        if (RESULT[J] >= THRESH) NFAIL++;
      }

      if (NFAIL > 0) NTESTF++;
      if (NTESTF == 1) {
        _printFirst(NOUNIT, PATH, THRESH);
        NTESTF = 2;
      }
      for (J = 1; J <= 17; J++) {
        if (RESULT[J] >= THRESH) {
          NOUNIT.println(
              ' N=${N.i5}, input example =${JTYPE.i3},  test(${J.i2})=${RESULT[J].g10_3}');
        }
      }

      NERRS += NFAIL;
      NTESTT += NTEST;
    }
  } catch (_) {}

  // Summary

  dlasum(PATH, NOUNIT, NERRS, NTESTT);
}

void _printFirst(Nout nout, String path, double threshold) {
  nout.println(
      '\n ${path.a3} -- Complex Schur Form Decomposition Expert Driver\n Matrix types (see ZDRVSX for details): ');
  nout.println(
      '\n Special Matrices:\n  1=Zero matrix.                          5=Diagonal: geometr. spaced entries.\n  2=Identity matrix.                      6=Diagonal: clustered entries.\n  3=Transposed Jordan block.              7=Diagonal: large, evenly spaced.\n  4=Diagonal: evenly spaced entries.      8=Diagonal: small, evenly spaced.');
  nout.println(
      ' Dense, Non-Symmetric Matrices:\n  9=Well-cond., evenly spaced eigenvals. 14=Ill-cond., geomet. spaced eigenals.\n 10=Well-cond., geom. spaced eigenvals.  15=Ill-conditioned, clustered e.vals.\n 11=Well-conditioned, clustered e.vals.  16=Ill-cond., random complex \n 12=Well-cond., random complex           17=Ill-cond., large rand. complx \n 13=Ill-conditioned, evenly spaced.      18=Ill-cond., small rand. complx ');
  nout.println(
      ' 19=Matrix with random O(1) entries.     21=Matrix with small random entries.\n 20=Matrix with large random entries.\n');
  nout.println(
      ' Tests performed with test threshold =${threshold.f8_2}\n ( A denotes A on input and T denotes A on output)\n\n 1 = 0 if T in Schur form (no sort),   1/ulp otherwise\n 2 = | A - VS T transpose(VS) | / ( n |A| ulp ) (no sort)\n 3 = | I - VS transpose(VS) | / ( n ulp ) (no sort) \n 4 = 0 if W are eigenvalues of T (no sort),  1/ulp otherwise\n 5 = 0 if T same no matter if VS computed (no sort),  1/ulp otherwise\n 6 = 0 if W same no matter if VS computed (no sort),  1/ulp otherwise');
  nout.println(
      ' 7 = 0 if T in Schur form (sort),   1/ulp otherwise\n 8 = | A - VS T transpose(VS) | / ( n |A| ulp ) (sort)\n 9 = | I - VS transpose(VS) | / ( n ulp ) (sort) \n 10 = 0 if W are eigenvalues of T (sort),  1/ulp otherwise\n 11 = 0 if T same no matter what else computed (sort),  1/ulp otherwise\n 12 = 0 if W same no matter what else computed (sort), 1/ulp otherwise\n 13 = 0 if sorting successful, 1/ulp otherwise\n 14 = 0 if RCONDE same no matter what else computed, 1/ulp otherwise\n 15 = 0 if RCONDv same no matter what else computed, 1/ulp otherwise\n 16 = | RCONDE - RCONDE(precomputed) | / cond(RCONDE),\n 17 = | RCONDV - RCONDV(precomputed) | / cond(RCONDV),');
}
