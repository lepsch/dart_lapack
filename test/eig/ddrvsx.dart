// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/range.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:test/test.dart';

import '../matgen/dlatme.dart';
import '../matgen/dlatmr.dart';
import '../matgen/dlatms.dart';
import '../test_driver.dart';
import 'common.dart';
import 'dget24.dart';
import 'dlasum.dart';

Future<void> ddrvsx(
  final int NSIZES,
  final Array<int> NN_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final Array<int> ISEED_,
  final double THRESH,
  final Nin NIUNIT,
  final Nout NOUNIT,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> H_,
  final Matrix<double> HT_,
  final Array<double> WR_,
  final Array<double> WI_,
  final Array<double> WRT_,
  final Array<double> WIT_,
  final Array<double> WRTMP_,
  final Array<double> WITMP_,
  final Matrix<double> VS_,
  final int LDVS,
  final Matrix<double> VS1_,
  final Array<double> RESULT_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final Array<bool> BWORK_,
  final Box<int> INFO,
  final TestDriver test,
  final String group,
) async {
  final NN = NN_.having();
  final DOTYPE = DOTYPE_.having();
  final ISEED = ISEED_.having();
  final A = A_.having(ld: LDA);
  final H = H_.having(ld: LDA);
  final HT = HT_.having(ld: LDA);
  final WR = WR_.having();
  final WI = WI_.having();
  final WRT = WRT_.having();
  final WIT = WIT_.having();
  final WRTMP = WRTMP_.having();
  final WITMP = WITMP_.having();
  final VS = VS_.having(ld: LDVS);
  final VS1 = VS1_.having(ld: LDVS);
  final RESULT = RESULT_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  final BWORK = BWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  const MAXTYP = 21;
  final ISLCT = Array<int>(20);
  const KTYPE = [
    1, 2, 3, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9, //
  ];
  const KMAGN = [
    1, 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 1, 2, 3, //
  ];
  const KMODE = [
    0, 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1, //
  ];
  const KCONDS = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0, 0, 0, //
  ];

  final PATH = '${'Double precision'[0]}SX';

  // Check for errors

  var NTESTT = 0;
  var NTESTF = 0;
  INFO.value = 0;

  // Important constants
  {
    var BADNN = false;

    // 12 is the largest dimension in the input file of precomputed
    // problems
    var NMAX = 12;
    for (var J = 1; J <= NSIZES; J++) {
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
      xerbla('DDRVSX', -INFO.value);
      return;
    }
  }

  var NERRS = 0;
  final PARAMS = claenv.IPARMS.copy();

  // If nothing to do check on NIUNIT

  if (NSIZES != 0 && NTYPES != 0) {
    // More Important constants

    final UNFL = dlamch('Safe minimum');
    final OVFL = ONE / UNFL;
    final ULP = dlamch('Precision');
    final ULPINV = ONE / ULP;
    final RTULP = sqrt(ULP);
    final RTULPI = ONE / RTULP;

    // Loop over sizes, types

    test.group(group, () {
      test.setUp(() {
        claenv.IPARMS.assign(PARAMS);
      });

      for (final JSIZE in 1.through(NSIZES)) {
        final N = NN[JSIZE];
        final MTYPES =
            NSIZES != 1 ? min(MAXTYP, NTYPES) : min(MAXTYP + 1, NTYPES);

        for (final JTYPE in 1.through(MTYPES)) {
          final skip = !DOTYPE[JTYPE];
          test('DDRVSX (N=$N, TYPE=$JTYPE)', () {
            // Save ISEED in case of an error.
            final IOLDSD = ISEED.copy();
            final IINFO = Box(0);

            // Compute "A"
            //
            // Control parameters:
            //
            // KMAGN  KCONDS  KMODE        KTYPE
            //    =1  O(1)   1       clustered 1  zero
            //    =2  large  large   clustered 2  identity
            //    =3  small          exponential  Jordan
            //    =4                 arithmetic   diagonal, (w/ eigenvalues)
            //    =5                 random log   symmetric, w/ eigenvalues
            //    =6                 random       general, w/ eigenvalues
            //    =7                              random diagonal
            //    =8                              random symmetric
            //    =9                              random general
            //    =10                             random triangular

            if (MTYPES <= MAXTYP) {
              final ITYPE = KTYPE[JTYPE - 1];
              final IMODE = KMODE[JTYPE - 1];

              // Compute norm
              final ANORM = switch (KMAGN[JTYPE - 1]) {
                1 => ONE,
                2 => OVFL * ULP,
                3 => UNFL * ULPINV,
                _ => throw UnimplementedError(),
              };

              dlaset('Full', LDA, N, ZERO, ZERO, A, LDA);
              IINFO.value = 0;
              final COND = ULPINV;

              // Special Matrices -- Identity & Jordan block

              // Zero

              if (ITYPE == 1) {
                IINFO.value = 0;
              } else if (ITYPE == 2) {
                // Identity

                for (var JCOL = 1; JCOL <= N; JCOL++) {
                  A[JCOL][JCOL] = ANORM;
                }
              } else if (ITYPE == 3) {
                // Jordan Block

                for (var JCOL = 1; JCOL <= N; JCOL++) {
                  A[JCOL][JCOL] = ANORM;
                  if (JCOL > 1) A[JCOL][JCOL - 1] = ONE;
                }
              } else if (ITYPE == 4) {
                // Diagonal Matrix, [Eigen]values Specified

                dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0,
                    'N', A, LDA, WORK(N + 1), IINFO);
              } else if (ITYPE == 5) {
                // Symmetric, eigenvalues specified

                dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, N, N,
                    'N', A, LDA, WORK(N + 1), IINFO);
              } else if (ITYPE == 6) {
                // General, eigenvalues specified

                final CONDS = switch (KCONDS[JTYPE - 1]) {
                  1 => ONE,
                  2 => RTULPI,
                  _ => ZERO,
                };

                final ADUMMA = ' ';
                dlatme(
                    N,
                    'S',
                    ISEED,
                    WORK,
                    IMODE,
                    COND,
                    ONE,
                    ADUMMA,
                    'T',
                    'T',
                    'T',
                    WORK(N + 1),
                    4,
                    CONDS,
                    N,
                    N,
                    ANORM,
                    A,
                    LDA,
                    WORK(2 * N + 1),
                    IINFO);
              } else if (ITYPE == 7) {
                // Diagonal, random eigenvalues

                final IDUMMA = Array<int>(1);
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
                    'NO',
                    A,
                    LDA,
                    IWORK,
                    IINFO);
              } else if (ITYPE == 8) {
                // Symmetric, random eigenvalues

                final IDUMMA = Array<int>(1);
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
                    N,
                    N,
                    ZERO,
                    ANORM,
                    'NO',
                    A,
                    LDA,
                    IWORK,
                    IINFO);
              } else if (ITYPE == 9) {
                // General, random eigenvalues

                final IDUMMA = Array<int>(1);
                dlatmr(
                    N,
                    N,
                    'S',
                    ISEED,
                    'N',
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
                    N,
                    N,
                    ZERO,
                    ANORM,
                    'NO',
                    A,
                    LDA,
                    IWORK,
                    IINFO);
                if (N >= 4) {
                  dlaset('Full', 2, N, ZERO, ZERO, A, LDA);
                  dlaset('Full', N - 3, 1, ZERO, ZERO, A(3, 1), LDA);
                  dlaset('Full', N - 3, 2, ZERO, ZERO, A(3, N - 1), LDA);
                  dlaset('Full', 1, N, ZERO, ZERO, A(N, 1), LDA);
                }
              } else if (ITYPE == 10) {
                // Triangular, random eigenvalues

                final IDUMMA = Array<int>(1);
                dlatmr(
                    N,
                    N,
                    'S',
                    ISEED,
                    'N',
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
                    N,
                    0,
                    ZERO,
                    ANORM,
                    'NO',
                    A,
                    LDA,
                    IWORK,
                    IINFO);
              } else {
                IINFO.value = 1;
              }

              if (IINFO.value != 0) {
                NOUNIT.println(
                    ' DDRVSX: Generator returned INFO=${IINFO.value.i6}.\n${' ' * 9}N=${N.i6}, JTYPE=${JTYPE.i6}, ISEED=(${IOLDSD.i5(4, ',')})');
                INFO.value = (IINFO.value).abs();
                return;
              }
            }

            // Test for minimal and generous workspace

            for (var IWK = 1; IWK <= 2; IWK++) {
              final NNWORK = max(IWK == 1 ? 3 * N : max(3 * N, 2 * N * N), 1);

              const NSLCT = 0;
              const RCDEIN = ZERO, RCDVIN = ZERO;
              dget24(
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
                  WR,
                  WI,
                  WRT,
                  WIT,
                  WRTMP,
                  WITMP,
                  VS,
                  LDVS,
                  VS1,
                  RCDEIN,
                  RCDVIN,
                  NSLCT,
                  ISLCT,
                  RESULT,
                  WORK,
                  NNWORK,
                  IWORK,
                  BWORK,
                  INFO);

              // Check for RESULT[j] > THRESH

              var NTEST = 0, NFAIL = 0;
              for (var J = 1; J <= 15; J++) {
                if (RESULT[J] >= ZERO) NTEST++;
                if (RESULT[J] >= THRESH) NFAIL++;
              }

              if (NFAIL > 0) NTESTF++;
              if (NTESTF == 1) {
                _printTestFailed(NOUNIT, PATH, THRESH);
                NTESTF = 2;
              }

              for (var J = 1; J <= 15; J++) {
                final reason =
                    ' N=${N.i5}, IWK=${IWK.i2}, seed=${IOLDSD.i4(4, ',')} type ${JTYPE.i2}, test(${J.i2})=${RESULT[J].g10_3}';
                test.expect(RESULT[J], lessThan(THRESH), reason: reason);
                if (RESULT[J] >= THRESH) {
                  NOUNIT.println(reason);
                }
              }

              NERRS += NFAIL;
              NTESTT += NTEST;
            }
          }, skip: skip);
        }
      }
    });
  }

  // Read in data from file to check accuracy of condition estimation
  // Read input data until N=0

  var JTYPE = 0;
  while (true) {
    final int N, NSLCT;
    final double RCDEIN, RCDVIN;
    try {
      (N, NSLCT) = await NIUNIT.readInt2();
      if (N == 0) break;
      JTYPE++;
      ISEED[1] = JTYPE;
      if (NSLCT > 0) await NIUNIT.readArray(ISLCT, NSLCT);
      await NIUNIT.readMatrix(A, N, N);
      (RCDEIN, RCDVIN) = await NIUNIT.readDouble2();
    } on EOF catch (_) {
      break;
    }

    final ctx = (
      A: A.copy(),
      ISEED: ISEED.copy(),
      JTYPE: JTYPE,
      ISLCT: ISLCT.copy(),
    );
    test.group(group, () {
      final (:A, :ISEED, :JTYPE, :ISLCT) = ctx;

      test.setUp(() {
        claenv.IPARMS.assign(PARAMS);
      });

      test('DDRVSX', () {
        dget24(
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
            WR,
            WI,
            WRT,
            WIT,
            WRTMP,
            WITMP,
            VS,
            LDVS,
            VS1,
            RCDEIN,
            RCDVIN,
            NSLCT,
            ISLCT,
            RESULT,
            WORK,
            LWORK,
            IWORK,
            BWORK,
            INFO);

        // Check for RESULT[j] > THRESH

        var NTEST = 0, NFAIL = 0;
        for (var J = 1; J <= 17; J++) {
          if (RESULT[J] >= ZERO) NTEST++;
          if (RESULT[J] >= THRESH) NFAIL++;
        }

        if (NFAIL > 0) NTESTF++;
        if (NTESTF == 1) {
          _printTestFailed(NOUNIT, PATH, THRESH);
          NTESTF = 2;
        }
        for (var J = 1; J <= 17; J++) {
          final reason =
              ' N=${N.i5}, input example =${JTYPE.i3},  test(${J.i2})=${RESULT[J].g10_3}';
          test.expect(RESULT[J], lessThan(THRESH), reason: reason);
          if (RESULT[J] >= THRESH) {
            NOUNIT.println(reason);
          }
        }

        NERRS += NFAIL;
        NTESTT += NTEST;
      });
    });
  }

  // Summary
  dlasum(PATH, NOUNIT, NERRS, NTESTT);
}

void _printTestFailed(
  final Nout nout,
  final String path,
  final double threshold,
) {
  nout.println(
      '\n ${path.a3} -- Real Schur Form Decomposition Expert Driver\n Matrix types (see DDRVSX for details):');
  nout.println(
      '\n Special Matrices:\n  1=Zero matrix.                          5=Diagonal: geometr. spaced entries.\n  2=Identity matrix.                      6=Diagonal: clustered entries.\n  3=Transposed Jordan block.              7=Diagonal: large, evenly spaced.\n  4=Diagonal: evenly spaced entries.      8=Diagonal: small, evenly spaced.');
  nout.println(
      ' Dense, Non-Symmetric Matrices:\n  9=Well-cond., evenly spaced eigenvals. 14=Ill-cond., geomet. spaced eigenals.\n 10=Well-cond., geom. spaced eigenvals.  15=Ill-conditioned, clustered e.vals.\n 11=Well-conditioned, clustered e.vals.  16=Ill-cond., random complex \n 12=Well-cond., random complex           17=Ill-cond., large rand. complx \n 13=Ill-conditioned, evenly spaced.      18=Ill-cond., small rand. complx ');
  nout.println(
      ' 19=Matrix with random O(1) entries.     21=Matrix with small random entries.\n 20=Matrix with large random entries.\n');
  nout.println(
      ' Tests performed with test threshold =${threshold.f8_2}\n ( A denotes A on input and T denotes A on output)\n\n 1 = 0 if T in Schur form (no sort),   1/ulp otherwise\n 2 = | A - VS T transpose(VS) | / ( n |A| ulp ) (no sort)\n 3 = | I - VS transpose(VS) | / ( n ulp ) (no sort) \n 4 = 0 if WR+sqrt(-1)*WI are eigenvalues of T (no sort),  1/ulp otherwise\n 5 = 0 if T same no matter if VS computed (no sort),  1/ulp otherwise\n 6 = 0 if WR, WI same no matter if VS computed (no sort),  1/ulp otherwise');
  nout.println(
      ' 7 = 0 if T in Schur form (sort),   1/ulp otherwise\n 8 = | A - VS T transpose(VS) | / ( n |A| ulp ) (sort)\n 9 = | I - VS transpose(VS) | / ( n ulp ) (sort) \n 10 = 0 if WR+sqrt(-1)*WI are eigenvalues of T (sort),  1/ulp otherwise\n 11 = 0 if T same no matter what else computed (sort),  1/ulp otherwise\n 12 = 0 if WR, WI same no matter what else computed (sort), 1/ulp otherwise\n 13 = 0 if sorting successful, 1/ulp otherwise\n 14 = 0 if RCONDE same no matter what else computed, 1/ulp otherwise\n 15 = 0 if RCONDv same no matter what else computed, 1/ulp otherwise\n 16 = | RCONDE - RCONDE(precomputed) | / cond(RCONDE),\n 17 = | RCONDV - RCONDV(precomputed) | / cond(RCONDV),');
}
