// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dnrm2.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgeev.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlapy2.dart';
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
import 'dget22.dart';
import 'dlasum.dart';

void ddrvev(
  final int NSIZES,
  final Array<int> NN_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final Array<int> ISEED_,
  final double THRESH,
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
  final Array<double> RESULT_,
  final Array<double> WORK_,
  final int NWORK,
  final Array<int> IWORK_,
  final Box<int> INFO,
  final TestDriver test,
) {
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
  final RESULT = RESULT_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  const TWO = 2.0;
  const MAXTYP = 21;
  final IDUMMA = Array<int>(1);
  final DUM = Array<double>(1), RES = Array<double>(2);
  final KTYPE = [
    1, 2, 3, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9 //
  ];
  final KMAGN = [
    1, 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 1, 2, 3 //
  ];
  final KMODE = [
    0, 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1 //
  ];
  final KCONDS = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0, 0, 0 //
  ];

  final PATH = '${'Double precision'[0]}EV';

  // Check for errors

  var NTESTT = 0, NTESTF = 0;
  INFO.value = 0;

  // Important constants
  {
    var BADNN = false;
    var NMAX = 0;
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
      // } else if ( NOUNIT <= 0 ) {
      //    INFO = -7;
    } else if (LDA < 1 || LDA < NMAX) {
      INFO.value = -9;
    } else if (LDVL < 1 || LDVL < NMAX) {
      INFO.value = -16;
    } else if (LDVR < 1 || LDVR < NMAX) {
      INFO.value = -18;
    } else if (LDLRE < 1 || LDLRE < NMAX) {
      INFO.value = -20;
    } else if (5 * NMAX + 2 * pow(NMAX, 2) > NWORK) {
      INFO.value = -23;
    }

    if (INFO.value != 0) {
      xerbla('DDRVEV', -INFO.value);
      return;
    }
  }

  // Quick return if nothing to do

  if (NSIZES == 0 || NTYPES == 0) return;

  // More Important constants

  final UNFL = dlamch('Safe minimum');
  final OVFL = ONE / UNFL;
  final ULP = dlamch('Precision');
  final ULPINV = ONE / ULP;
  final RTULP = sqrt(ULP);
  final RTULPI = ONE / RTULP;

  // Loop over sizes, types

  var NERRS = 0;

  for (final JSIZE in 1.through(NSIZES)) {
    final N = NN[JSIZE];
    final MTYPES = NSIZES != 1 ? min(MAXTYP, NTYPES) : min(MAXTYP + 1, NTYPES);

    for (final JTYPE in 1.through(MTYPES)) {
      final skip = !DOTYPE[JTYPE];
      test('DDRVEV (N=$N, TYPE=$JTYPE)', () {
        // Save ISEED in case of an error.
        final IOLDSD = ISEED.copy();
        final IINFO = Box(0);

        // Compute "A"

        // Control parameters:

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

            dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'N',
                A, LDA, WORK(N + 1), IINFO);
          } else if (ITYPE == 5) {
            // Symmetric, eigenvalues specified

            dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, N, N, 'N',
                A, LDA, WORK(N + 1), IINFO);
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

          test.expect(IINFO.value, 0);
          if (IINFO.value != 0) {
            _print9993(NOUNIT, 'Generator', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            return;
          }
        }

        // Test for minimal and generous workspace

        for (var IWK = 1; IWK <= 2; IWK++) {
          final NNWORK =
              max(IWK == 1 ? 4 * N : 5 * N + 2 * pow(N, 2).toInt(), 1);

          // Initialize RESULT

          for (var J = 1; J <= 7; J++) {
            RESULT[J] = -ONE;
          }

          // Compute eigenvalues and eigenvectors, and test them

          dlacpy('F', N, N, A, LDA, H, LDA);
          dgeev('V', 'V', N, H, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, NNWORK,
              IINFO);
          if (IINFO.value != 0) {
            RESULT[1] = ULPINV;
            _print9993(NOUNIT, 'DGEEV1', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            break;
          }

          // Do Test (1)

          dget22('N', 'N', 'N', N, A, LDA, VR, LDVR, WR, WI, WORK, RES);
          RESULT[1] = RES[1];

          // Do Test (2)

          dget22('T', 'N', 'T', N, A, LDA, VL, LDVL, WR, WI, WORK, RES);
          RESULT[2] = RES[1];

          // Do Test (3)

          for (var J = 1; J <= N; J++) {
            final double TNRM;
            if (WI[J] == ZERO) {
              TNRM = dnrm2(N, VR(1, J).asArray(), 1);
            } else if (WI[J] > ZERO) {
              TNRM = dlapy2(dnrm2(N, VR(1, J).asArray(), 1),
                  dnrm2(N, VR(1, J + 1).asArray(), 1));
            } else {
              TNRM = ONE;
            }
            RESULT[3] = max(RESULT[3], min(ULPINV, (TNRM - ONE).abs() / ULP));
            if (WI[J] > ZERO) {
              var VMX = ZERO;
              var VRMX = ZERO;
              for (var JJ = 1; JJ <= N; JJ++) {
                final VTST = dlapy2(VR[JJ][J], VR[JJ][J + 1]);
                if (VTST > VMX) VMX = VTST;
                if (VR[JJ][J + 1] == ZERO && VR[JJ][J].abs() > VRMX) {
                  VRMX = VR[JJ][J].abs();
                }
              }
              if (VRMX / VMX < ONE - TWO * ULP) RESULT[3] = ULPINV;
            }
          }

          // Do Test (4)

          for (var J = 1; J <= N; J++) {
            final double TNRM;
            if (WI[J] == ZERO) {
              TNRM = dnrm2(N, VL(1, J).asArray(), 1);
            } else if (WI[J] > ZERO) {
              TNRM = dlapy2(dnrm2(N, VL(1, J).asArray(), 1),
                  dnrm2(N, VL(1, J + 1).asArray(), 1));
            } else {
              TNRM = ONE;
            }
            RESULT[4] = max(RESULT[4], min(ULPINV, (TNRM - ONE).abs() / ULP));
            if (WI[J] > ZERO) {
              var VMX = ZERO;
              var VRMX = ZERO;
              for (var JJ = 1; JJ <= N; JJ++) {
                final VTST = dlapy2(VL[JJ][J], VL[JJ][J + 1]);
                if (VTST > VMX) VMX = VTST;
                if (VL[JJ][J + 1] == ZERO && VL[JJ][J].abs() > VRMX) {
                  VRMX = VL[JJ][J].abs();
                }
              }
              if (VRMX / VMX < ONE - TWO * ULP) RESULT[4] = ULPINV;
            }
          }

          // Compute eigenvalues only, and test them
          while (true) {
            dlacpy('F', N, N, A, LDA, H, LDA);
            dgeev('N', 'N', N, H, LDA, WR1, WI1, DUM.asMatrix(1), 1,
                DUM.asMatrix(1), 1, WORK, NNWORK, IINFO);
            if (IINFO.value != 0) {
              RESULT[1] = ULPINV;
              _print9993(NOUNIT, 'DGEEV2', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              break;
            }

            // Do Test (5)

            for (var J = 1; J <= N; J++) {
              if (WR[J] != WR1[J] || WI[J] != WI1[J]) RESULT[5] = ULPINV;
            }

            // Compute eigenvalues and right eigenvectors, and test them

            dlacpy('F', N, N, A, LDA, H, LDA);
            dgeev('N', 'V', N, H, LDA, WR1, WI1, DUM.asMatrix(1), 1, LRE, LDLRE,
                WORK, NNWORK, IINFO);
            if (IINFO.value != 0) {
              RESULT[1] = ULPINV;
              _print9993(NOUNIT, 'DGEEV3', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              break;
            }

            // Do Test (5) again

            for (var J = 1; J <= N; J++) {
              if (WR[J] != WR1[J] || WI[J] != WI1[J]) RESULT[5] = ULPINV;
            }

            // Do Test (6)

            for (var J = 1; J <= N; J++) {
              for (var JJ = 1; JJ <= N; JJ++) {
                if (VR[J][JJ] != LRE[J][JJ]) RESULT[6] = ULPINV;
              }
            }

            // Compute eigenvalues and left eigenvectors, and test them

            dlacpy('F', N, N, A, LDA, H, LDA);
            dgeev('V', 'N', N, H, LDA, WR1, WI1, LRE, LDLRE, DUM.asMatrix(1), 1,
                WORK, NNWORK, IINFO);
            if (IINFO.value != 0) {
              RESULT[1] = ULPINV;
              _print9993(NOUNIT, 'DGEEV4', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              break;
            }

            // Do Test (5) again

            for (var J = 1; J <= N; J++) {
              if (WR[J] != WR1[J] || WI[J] != WI1[J]) RESULT[5] = ULPINV;
            }

            // Do Test (7)

            for (var J = 1; J <= N; J++) {
              for (var JJ = 1; JJ <= N; JJ++) {
                if (VL[J][JJ] != LRE[J][JJ]) RESULT[7] = ULPINV;
              }
            }

            break;
          }

          // End of Loop -- Check for RESULT[j] > THRESH

          var NTEST = 0, NFAIL = 0;
          for (var J = 1; J <= 7; J++) {
            if (RESULT[J] >= ZERO) NTEST++;
            if (RESULT[J] >= THRESH) NFAIL++;
          }

          if (NFAIL > 0) NTESTF++;
          if (NTESTF == 1) {
            NOUNIT.println(
                '\n ${PATH.a3} -- Real Eigenvalue-Eigenvector Decomposition Driver\n Matrix types (see DDRVEV for details): ');
            NOUNIT.println(
                '\n Special Matrices:\n  1=Zero matrix.                          5=Diagonal: geometr. spaced entries.\n  2=Identity matrix.                      6=Diagonal: clustered entries.\n  3=Transposed Jordan block.              7=Diagonal: large, evenly spaced.\n  4=Diagonal: evenly spaced entries.      8=Diagonal: small, evenly spaced.');
            NOUNIT.println(
                ' Dense, Non-Symmetric Matrices:\n  9=Well-cond., evenly spaced eigenvals. 14=Ill-cond., geomet. spaced eigenals.\n 10=Well-cond., geom. spaced eigenvals.  15=Ill-conditioned, clustered e.vals.\n 11=Well-conditioned, clustered e.vals.  16=Ill-cond., random complex \n 12=Well-cond., random complex ${' ' * 6}    17=Ill-cond., large rand. complx \n 13=Ill-conditioned, evenly spaced.      18=Ill-cond., small rand. complx ');
            NOUNIT.println(
                ' 19=Matrix with random O(1) entries.     21=Matrix with small random entries.\n 20=Matrix with large random entries.\n');
            NOUNIT.println(
                ' Tests performed with test threshold =${THRESH.f8_2}\n\n 1 = | A VR - VR W | / ( n |A| ulp ) \n 2 = | transpose(A) VL - VL W | / ( n |A| ulp ) \n 3 = | |VR(i)| - 1 | / ulp \n 4 = | |VL(i)| - 1 | / ulp \n 5 = 0 if W same no matter if VR or VL computed, 1/ulp otherwise\n 6 = 0 if VR same no matter if VL computed,  1/ulp otherwise\n 7 = 0 if VL same no matter if VR computed,  1/ulp otherwise\n');
            NTESTF = 2;
          }

          for (var J = 1; J <= 7; J++) {
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

  // Summary
  dlasum(PATH, NOUNIT, NERRS, NTESTT);
}

void _print9993(
  final Nout nout,
  final String s,
  final int info,
  final int n,
  final int jtype,
  final Array<int> iseed,
) {
  nout.println(
      ' DDRVEV: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}
