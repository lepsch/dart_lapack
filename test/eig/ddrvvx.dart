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
import 'dget23.dart';
import 'dlasum.dart';

Future<void> ddrvvx(
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
  final Array<double> WR_,
  final Array<double> WI_,
  final Array<double> WR1_,
  final Array<double> WI1_,
  final Matrix<double> VL_,
  final int LDVL,
  final Matrix<double> VR_,
  final int LDVR,
  final Matrix<double> LRE_,
  final int LDLRE,
  final Array<double> RCONDV_,
  final Array<double> RCNDV1_,
  final Array<double> RCDVIN_,
  final Array<double> RCONDE_,
  final Array<double> RCNDE1_,
  final Array<double> RCDEIN_,
  final Array<double> SCALE_,
  final Array<double> SCALE1_,
  final Array<double> RESULT_,
  final Array<double> WORK_,
  final int NWORK,
  final Array<int> IWORK_,
  final Box<int> INFO,
  final TestDriver test,
  final String group,
) async {
  final NN = NN_.having();
  final DOTYPE = DOTYPE_.having();
  final ISEED = ISEED_.having();
  final A = A_.having(ld: LDA);
  final H = H_.having(ld: LDA);
  final WR = WR_.having();
  final WI = WI_.having();
  final WR1 = WR1_.having();
  final WI1 = WI1_.having();
  final VL = VL_.having(ld: LDVL);
  final VR = VR_.having(ld: LDVR);
  final LRE = LRE_.having(ld: LDLRE);
  final RCONDV = RCONDV_.having();
  final RCNDV1 = RCNDV1_.having();
  final RCDVIN = RCDVIN_.having();
  final RCONDE = RCONDE_.having();
  final RCNDE1 = RCNDE1_.having();
  final RCDEIN = RCDEIN_.having();
  final SCALE = SCALE_.having();
  final SCALE1 = SCALE1_.having();
  final RESULT = RESULT_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  const MAXTYP = 21;
  const KTYPE = [
    1, 2, 3, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9 //
  ];
  const KMAGN = [
    1, 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 1, 2, 3 //
  ];
  const KMODE = [
    0, 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1 //
  ];
  const KCONDS = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0, 0, 0 //
  ];
  const BAL = ['N', 'P', 'S', 'B'];

  final PATH = '${'Double precision'[0]}VX';

  // Check for errors

  var NTESTT = 0, NTESTF = 0;
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
    } else if (LDA < 1 || LDA < NMAX) {
      INFO.value = -10;
    } else if (LDVL < 1 || LDVL < NMAX) {
      INFO.value = -17;
    } else if (LDVR < 1 || LDVR < NMAX) {
      INFO.value = -19;
    } else if (LDLRE < 1 || LDLRE < NMAX) {
      INFO.value = -21;
    } else if (6 * NMAX + 2 * pow(NMAX, 2) > NWORK) {
      INFO.value = -32;
    }

    if (INFO.value != 0) {
      xerbla('DDRVVX', -INFO.value);
      return;
    }
  }

  final PARAMS = claenv.IPARMS.copy();

  var NERRS = 0;

  // If nothing to do check on NIUNIT

  if (NSIZES != 0 && NTYPES != 0) {
    // More Important constants
    final UNFL = dlamch('Safe minimum');
    final OVFL = ONE / UNFL;
    final ULP = dlamch('Precision');
    final ULPINV = ONE / ULP;
    final RTULP = sqrt(ULP);
    final RTULPI = ONE / RTULP;

    test.group(group, () {
      test.setUp(() {
        claenv.IPARMS.assign(PARAMS);
      });

      // Loop over sizes, types
      for (final JSIZE in 1.through(NSIZES)) {
        final N = NN[JSIZE];
        final MTYPES =
            NSIZES != 1 ? min(MAXTYP, NTYPES) : min(MAXTYP + 1, NTYPES);

        for (final JTYPE in 1.through(MTYPES)) {
          final skip = !DOTYPE[JTYPE];
          test('DDRVVX (N=$N, TYPE=$JTYPE)', () {
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
                    ' DDRVVX: Generator returned INFO=${IINFO.value.i6}.\n${' ' * 9}N=${N.i6}, JTYPE=${JTYPE.i6}, ISEED=(${IOLDSD.i5(4, ',')})');
                INFO.value = (IINFO.value).abs();
                return;
              }
            }

            // Test for minimal and generous workspace

            for (var IWK = 1; IWK <= 3; IWK++) {
              final NNWORK = max(
                  switch (IWK) {
                    1 => 3 * N,
                    2 => 6 * N + pow(N, 2).toInt(),
                    _ => 6 * N + 2 * pow(N, 2).toInt(),
                  },
                  1);

              // Test for all balancing options

              for (var IBAL = 1; IBAL <= 4; IBAL++) {
                final BALANC = BAL[IBAL - 1];

                // Perform tests

                dget23(
                    false,
                    BALANC,
                    JTYPE,
                    THRESH,
                    IOLDSD,
                    NOUNIT,
                    N,
                    A,
                    LDA,
                    H,
                    WR,
                    WI,
                    WR1,
                    WI1,
                    VL,
                    LDVL,
                    VR,
                    LDVR,
                    LRE,
                    LDLRE,
                    RCONDV,
                    RCNDV1,
                    RCDVIN,
                    RCONDE,
                    RCNDE1,
                    RCDEIN,
                    SCALE,
                    SCALE1,
                    RESULT,
                    WORK,
                    NNWORK,
                    IWORK,
                    INFO);

                // Check for RESULT[j] > THRESH

                var NTEST = 0, NFAIL = 0;
                for (var J = 1; J <= 9; J++) {
                  if (RESULT[J] >= ZERO) NTEST++;
                  if (RESULT[J] >= THRESH) NFAIL++;
                }

                if (NFAIL > 0) NTESTF++;
                if (NTESTF == 1) {
                  _printTestFailed(NOUNIT, PATH, THRESH);
                  NTESTF = 2;
                }

                for (var J = 1; J <= 9; J++) {
                  final reason =
                      ' BALANC=\'${BALANC.a1}\',N=${N.i4},IWK=${IWK.i1}, seed=${IOLDSD.i4(4, ',')} type ${JTYPE.i2}, test(${J.i2})=${RESULT[J].g10_3}';
                  test.expect(RESULT[J], lessThan(THRESH), reason: reason);
                  if (RESULT[J] >= THRESH) {
                    NOUNIT.println(reason);
                  }
                }

                NERRS += NFAIL;
                NTESTT += NTEST;
              }
            }
          }, skip: skip);
        }
      }
    });
  }

  // Read in data from file to check accuracy of condition estimation.
  // Assume input eigenvalues are sorted lexicographically (increasing
  // by real part, then decreasing by imaginary part)

  var JTYPE = 0;
  while (true) {
    final int N;
    try {
      N = await NIUNIT.readInt();

      // Read input data until N=0
      if (N == 0) break;

      JTYPE++;
      ISEED[1] = JTYPE;
      await NIUNIT.readMatrix(A, N, N);
      for (var I = 1; I <= N; I++) {
        final (f1, f2, f3, f4) = await NIUNIT.readDouble4();
        WR1[I] = f1;
        WI1[I] = f2;
        RCDEIN[I] = f3;
        RCDVIN[I] = f4;
      }
    } on EOF catch (_) {
      break;
    }

    final ctx = (
      A: A.copy(),
      WR1: WR1.copy(),
      WI1: WI1.copy(),
      RCDEIN: RCDEIN.copy(),
      RCDVIN: RCDVIN.copy(),
      ISEED: ISEED,
      JTYPE: JTYPE,
    );
    test.group(group, () {
      final (:A, :WR1, :WI1, :RCDEIN, :RCDVIN, :ISEED, :JTYPE) = ctx;

      test.setUp(() {
        claenv.IPARMS.assign(PARAMS);
      });

      test('DDRVVX - precomputed', () {
        dget23(
            true,
            'N',
            22,
            THRESH,
            ISEED,
            NOUNIT,
            N,
            A,
            LDA,
            H,
            WR,
            WI,
            WR1,
            WI1,
            VL,
            LDVL,
            VR,
            LDVR,
            LRE,
            LDLRE,
            RCONDV,
            RCNDV1,
            RCDVIN,
            RCONDE,
            RCNDE1,
            RCDEIN,
            SCALE,
            SCALE1,
            RESULT,
            WORK,
            6 * N + 2 * pow(N, 2).toInt(),
            IWORK,
            INFO);

        // Check for RESULT[j] > THRESH

        var NTEST = 0, NFAIL = 0;
        for (var J = 1; J <= 11; J++) {
          if (RESULT[J] >= ZERO) NTEST++;
          if (RESULT[J] >= THRESH) NFAIL++;
        }

        if (NFAIL > 0) NTESTF++;
        if (NTESTF == 1) {
          _printTestFailed(NOUNIT, PATH, THRESH);
          NTESTF = 2;
        }

        for (var J = 1; J <= 11; J++) {
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
      '\n ${path.a3} -- Real Eigenvalue-Eigenvector Decomposition Expert Driver\n Matrix types (see DDRVVX for details): ');
  nout.println(
      '\n Special Matrices:\n  1=Zero matrix.                          5=Diagonal: geometr. spaced entries.\n  2=Identity matrix.                      6=Diagonal: clustered entries.\n  3=Transposed Jordan block.              7=Diagonal: large, evenly spaced.\n  4=Diagonal: evenly spaced entries.      8=Diagonal: small, evenly spaced.');
  nout.println(
      ' Dense, Non-Symmetric Matrices:\n  9=Well-cond., evenly spaced eigenvals. 14=Ill-cond., geomet. spaced eigenals.\n 10=Well-cond., geom. spaced eigenvals.  15=Ill-conditioned, clustered e.vals.\n 11=Well-conditioned, clustered e.vals.  16=Ill-cond., random complex \n 12=Well-cond., random complex           17=Ill-cond., large rand. complx \n 13=Ill-conditioned, evenly spaced.      18=Ill-cond., small rand. complx ');
  nout.println(
      ' 19=Matrix with random O(1) entries.     21=Matrix with small random entries.\n 20=Matrix with large random entries.    22=Matrix read from input file\n');
  nout.println(
      ' Tests performed with test threshold =${threshold.f8_2}\n\n 1 = | A VR - VR W | / ( n |A| ulp ) \n 2 = | transpose(A) VL - VL W | / ( n |A| ulp ) \n 3 = | |VR(i)| - 1 | / ulp \n 4 = | |VL(i)| - 1 | / ulp \n 5 = 0 if W same no matter if VR or VL computed, 1/ulp otherwise\n 6 = 0 if VR same no matter what else computed,  1/ulp otherwise\n 7 = 0 if VL same no matter what else computed,  1/ulp otherwise\n 8 = 0 if RCONDV same no matter what else computed,  1/ulp otherwise\n 9 = 0 if SCALE, ILO, IHI, ABNRM same no matter what else computed,  1/ulp otherwise\n 10 = | RCONDV - RCONDV[precomputed] | / cond(RCONDV),\n 11 = | RCONDE - RCONDE[precomputed] | / cond(RCONDE),');
}
