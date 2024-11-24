// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgehrd.dart';
import 'package:lapack/src/dhsein.dart';
import 'package:lapack/src/dhseqr.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dorghr.dart';
import 'package:lapack/src/dormhr.dart';
import 'package:lapack/src/dtrevc.dart';
import 'package:lapack/src/dtrevc3.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/range.dart';
import 'package:lapack/src/xerbla.dart';

import '../matgen/dlatme.dart';
import '../matgen/dlatmr.dart';
import '../matgen/dlatms.dart';
import '../test_driver.dart';
import 'dlafts.dart';
import 'dget10.dart';
import 'dget22.dart';
import 'dhst01.dart';
import 'dlasum.dart';

void dchkhs(
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
  final Matrix<double> T1_,
  final Matrix<double> T2_,
  final Matrix<double> U_,
  final int LDU,
  final Matrix<double> Z_,
  final Matrix<double> UZ_,
  final Array<double> WR1_,
  final Array<double> WI1_,
  final Array<double> WR2_,
  final Array<double> WI2_,
  final Array<double> WR3_,
  final Array<double> WI3_,
  final Matrix<double> EVECTL_,
  final Matrix<double> EVECTR_,
  final Matrix<double> EVECTY_,
  final Matrix<double> EVECTX_,
  final Matrix<double> UU_,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final int NWORK,
  final Array<int> IWORK_,
  final Array<bool> SELECT_,
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
  final H = H_.having(ld: LDA);
  final T1 = T1_.having(ld: LDA);
  final T2 = T2_.having(ld: LDA);
  final U = U_.having(ld: LDU);
  final Z = Z_.having(ld: LDU);
  final UZ = UZ_.having(ld: LDU);
  final WR1 = WR1_.having();
  final WI1 = WI1_.having();
  final WR2 = WR2_.having();
  final WI2 = WI2_.having();
  final WR3 = WR3_.having();
  final WI3 = WI3_.having();
  final EVECTL = EVECTL_.having(ld: LDU);
  final EVECTR = EVECTR_.having(ld: LDU);
  final EVECTY = EVECTY_.having(ld: LDU);
  final EVECTX = EVECTX_.having(ld: LDU);
  final UU = UU_.having(ld: LDU);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  final SELECT = SELECT_.having();
  final RESULT = RESULT_.having();
  const ZERO = 0.0, ONE = 1.0;
  const MAXTYP = 21;
  final IDUMMA = Array<int>(1);
  final DUMMA = Array<double>(6);
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

  // Check for errors

  var NTESTT = 0;
  INFO.value = 0;

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
  } else if (LDA <= 1 || LDA < NMAX) {
    INFO.value = -9;
  } else if (LDU <= 1 || LDU < NMAX) {
    INFO.value = -14;
  } else if (4 * NMAX * NMAX + 2 > NWORK) {
    INFO.value = -28;
  }

  if (INFO.value != 0) {
    xerbla('DCHKHS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (NSIZES == 0 || NTYPES == 0) return;

  // More important constants

  final UNFL = dlamch('Safe minimum');
  final OVFL = dlamch('Overflow');
  final ULP = dlamch('Epsilon') * dlamch('Base');
  final ULPINV = ONE / ULP;
  final RTUNFL = sqrt(UNFL);
  final RTOVFL = sqrt(OVFL);
  final RTULP = sqrt(ULP);
  final RTULPI = ONE / RTULP;

  // Loop over sizes, types

  final NERRS = Box(0);

  for (final JSIZE in 1.through(NSIZES)) {
    final N = NN[JSIZE];
    if (N == 0) continue;
    final N1 = max(1, N);
    final ANINV = ONE / N1;
    final MTYPES = NSIZES != 1 ? min(MAXTYP, NTYPES) : min(MAXTYP + 1, NTYPES);

    for (final JTYPE in 1.through(MTYPES)) {
      final skip = !DOTYPE[JTYPE];
      test('DCHKHS (SIZE = $N, TYPE = $JTYPE)', () {
        var NTEST = 0;
        final IINFO = Box(0), IN = Box(0);

        // Save ISEED in case of an error.

        final IOLDSD = ISEED.copy();

        // Initialize RESULT

        for (var J = 1; J <= 16; J++) {
          RESULT[J] = ZERO;
        }

        // Compute "A"
        //
        // Control parameters:
        //
        //     KMAGN  KCONDS  KMODE        KTYPE
        //        =1  O(1)   1       clustered 1  zero
        //        =2  large  large   clustered 2  identity
        //        =3  small          exponential  Jordan
        //        =4                 arithmetic   diagonal, (w/ eigenvalues)
        //        =5                 random log   symmetric, w/ eigenvalues
        //        =6                 random       general, w/ eigenvalues
        //        =7                              random diagonal
        //        =8                              random symmetric
        //        =9                              random general
        //        =10                             random triangular

        // Compute norm
        final ANORM = switch (KMAGN[JTYPE - 1]) {
          1 => ONE,
          2 => (RTOVFL * ULP) * ANINV,
          3 => RTUNFL * N * ULPINV,
          _ => throw UnimplementedError(),
        };

        if (MTYPES <= MAXTYP) {
          final ITYPE = KTYPE[JTYPE - 1];
          final IMODE = KMODE[JTYPE - 1];

          dlaset('Full', LDA, N, ZERO, ZERO, A, LDA);
          IINFO.value = 0;
          final COND = ULPINV;

          // Special Matrices

          if (ITYPE == 1) {
            // Zero

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

            const ADUMMA = ' ';
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
            print9999(NOUNIT, 'Generator', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            return;
          }
        }

        while (true) {
          // Call DGEHRD to compute H and U, do tests.

          dlacpy(' ', N, N, A, LDA, H, LDA);

          NTEST = 1;

          final ILO = 1, IHI = N;

          dgehrd(N, ILO, IHI, H, LDA, WORK, WORK(N + 1), NWORK - N, IINFO);
          if (IINFO.value != 0) {
            RESULT[1] = ULPINV;
            print9999(NOUNIT, 'DGEHRD', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            break;
          }

          for (var J = 1; J <= N - 1; J++) {
            UU[J + 1][J] = ZERO;
            for (var I = J + 2; I <= N; I++) {
              U[I][J] = H[I][J];
              UU[I][J] = H[I][J];
              H[I][J] = ZERO;
            }
          }
          dcopy(N - 1, WORK, 1, TAU, 1);
          dorghr(N, ILO, IHI, U, LDU, WORK, WORK(N + 1), NWORK - N, IINFO);
          NTEST = 2;

          dhst01(N, ILO, IHI, A, LDA, H, LDA, U, LDU, WORK, NWORK, RESULT(1));

          // Call DHSEQR to compute T1, T2 and Z, do tests.

          // Eigenvalues only (WR3,WI3)

          dlacpy(' ', N, N, H, LDA, T2, LDA);
          NTEST = 3;
          RESULT[3] = ULPINV;

          dhseqr('E', 'N', N, ILO, IHI, T2, LDA, WR3, WI3, UZ, LDU, WORK, NWORK,
              IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DHSEQR(E)', IINFO.value, N, JTYPE, IOLDSD);
            if (IINFO.value <= N + 2) {
              INFO.value = IINFO.value.abs();
              break;
            }
          }

          // Eigenvalues (WR2,WI2) and Full Schur Form (T2)

          dlacpy(' ', N, N, H, LDA, T2, LDA);

          dhseqr('S', 'N', N, ILO, IHI, T2, LDA, WR2, WI2, UZ, LDU, WORK, NWORK,
              IINFO);
          if (IINFO.value != 0 && IINFO.value <= N + 2) {
            print9999(NOUNIT, 'DHSEQR(S)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            break;
          }

          // Eigenvalues (WR1,WI1), Schur Form (T1), and Schur vectors
          // (UZ)

          dlacpy(' ', N, N, H, LDA, T1, LDA);
          dlacpy(' ', N, N, U, LDU, UZ, LDU);

          dhseqr('S', 'V', N, ILO, IHI, T1, LDA, WR1, WI1, UZ, LDU, WORK, NWORK,
              IINFO);
          if (IINFO.value != 0 && IINFO.value <= N + 2) {
            print9999(NOUNIT, 'DHSEQR(V)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            break;
          }

          // Compute Z = U' UZ

          dgemm('T', 'N', N, N, N, ONE, U, LDU, UZ, LDU, ZERO, Z, LDU);
          NTEST = 8;

          // Do Tests 3: | H - Z T Z' | / ( |H| n ulp )
          // and 4: | I - Z Z' | / ( n ulp )

          dhst01(N, ILO, IHI, H, LDA, T1, LDA, Z, LDU, WORK, NWORK, RESULT(3));

          // Do Tests 5: | A - UZ T (UZ)' | / ( |A| n ulp )
          // and 6: | I - UZ (UZ)' | / ( n ulp )

          dhst01(N, ILO, IHI, A, LDA, T1, LDA, UZ, LDU, WORK, NWORK, RESULT(5));

          // Do Test 7: | T2 - T1 | / ( |T| n ulp )

          dget10(N, N, T2, LDA, T1, LDA, WORK, RESULT.box(7));

          // Do Test 8: | W2 - W1 | / ( max(|W1|,|W2|) ulp )

          var TEMP1 = ZERO, TEMP2 = ZERO;
          for (var J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1,
                max(WR1[J].abs() + WI1[J].abs(), WR2[J].abs() + WI2[J].abs()));
            TEMP2 =
                max(TEMP2, (WR1[J] - WR2[J]).abs() + (WI1[J] - WI2[J]).abs());
          }

          RESULT[8] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          // Compute the Left and Right Eigenvectors of T

          // Compute the Right eigenvector Matrix:

          NTEST = 9;
          RESULT[9] = ULPINV;

          // Select last max(N/4,1) real, max(N/4,1) complex eigenvectors

          var NSELC = 0, NSELR = 0, J = N;
          do {
            if (WI1[J] == ZERO) {
              if (NSELR < max(N ~/ 4, 1)) {
                NSELR++;
                SELECT[J] = true;
              } else {
                SELECT[J] = false;
              }
              J--;
            } else {
              if (NSELC < max(N ~/ 4, 1)) {
                NSELC++;
                SELECT[J] = true;
                SELECT[J - 1] = false;
              } else {
                SELECT[J] = false;
                SELECT[J - 1] = false;
              }
              J -= 2;
            }
          } while (J > 0);

          dtrevc('Right', 'All', SELECT, N, T1, LDA, DUMMA.asMatrix(LDU), LDU,
              EVECTR, LDU, N, IN, WORK, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DTREVC(R,A)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            break;
          }

          // Test 9:  | TR - RW | / ( |T| |R| ulp )

          dget22('N', 'N', 'N', N, T1, LDA, EVECTR, LDU, WR1, WI1, WORK, DUMMA);
          RESULT[9] = DUMMA[1];
          if (DUMMA[2] > THRESH) {
            print9998(NOUNIT, 'Right', 'DTREVC', DUMMA[2], N, JTYPE, IOLDSD);
          }

          // Compute selected right eigenvectors and confirm that
          // they agree with previous right eigenvectors

          dtrevc('Right', 'Some', SELECT, N, T1, LDA, DUMMA.asMatrix(LDU), LDU,
              EVECTL, LDU, N, IN, WORK, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DTREVC(R,S)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            break;
          }

          var MATCH = true;
          matchRightLoop:
          for (var J = 1, K = 1; J <= N; J++) {
            if (SELECT[J] && WI1[J] == ZERO) {
              for (var JJ = 1; JJ <= N; JJ++) {
                if (EVECTR[JJ][J] != EVECTL[JJ][K]) {
                  MATCH = false;
                  break matchRightLoop;
                }
              }
              K++;
            } else if (SELECT[J] && WI1[J] != ZERO) {
              for (var JJ = 1; JJ <= N; JJ++) {
                if (EVECTR[JJ][J] != EVECTL[JJ][K] ||
                    EVECTR[JJ][J + 1] != EVECTL[JJ][K + 1]) {
                  MATCH = false;
                  break matchRightLoop;
                }
              }
              K += 2;
            }
          }

          if (!MATCH) print9997(NOUNIT, 'Right', 'DTREVC', N, JTYPE, IOLDSD);

          // Compute the Left eigenvector Matrix:

          NTEST = 10;
          RESULT[10] = ULPINV;
          dtrevc('Left', 'All', SELECT, N, T1, LDA, EVECTL, LDU,
              DUMMA.asMatrix(LDU), LDU, N, IN, WORK, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DTREVC(L,A)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            break;
          }

          // Test 10:  | LT - WL | / ( |T| |L| ulp )

          dget22('Trans', 'N', 'Conj', N, T1, LDA, EVECTL, LDU, WR1, WI1, WORK,
              DUMMA(3));
          RESULT[10] = DUMMA[3];
          if (DUMMA[4] > THRESH) {
            print9998(NOUNIT, 'Left', 'DTREVC', DUMMA[4], N, JTYPE, IOLDSD);
          }

          // Compute selected left eigenvectors and confirm that
          // they agree with previous left eigenvectors

          dtrevc('Left', 'Some', SELECT, N, T1, LDA, EVECTR, LDU,
              DUMMA.asMatrix(LDU), LDU, N, IN, WORK, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DTREVC(L,S)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            break;
          }

          MATCH = true;
          matchLeftLoop:
          for (var J = 1, K = 1; J <= N; J++) {
            if (SELECT[J] && WI1[J] == ZERO) {
              for (var JJ = 1; JJ <= N; JJ++) {
                if (EVECTL[JJ][J] != EVECTR[JJ][K]) {
                  MATCH = false;
                  break matchLeftLoop;
                }
              }
              K++;
            } else if (SELECT[J] && WI1[J] != ZERO) {
              for (var JJ = 1; JJ <= N; JJ++) {
                if (EVECTL[JJ][J] != EVECTR[JJ][K] ||
                    EVECTL[JJ][J + 1] != EVECTR[JJ][K + 1]) {
                  MATCH = false;
                  break matchLeftLoop;
                }
              }
              K += 2;
            }
          }

          if (!MATCH) print9997(NOUNIT, 'Left', 'DTREVC', N, JTYPE, IOLDSD);

          // Call DHSEIN for Right eigenvectors of H, do test 11

          NTEST = 11;
          RESULT[11] = ULPINV;
          for (var J = 1; J <= N; J++) {
            SELECT[J] = true;
          }

          dhsein(
              'Right',
              'Qr',
              'Ninitv',
              SELECT,
              N,
              H,
              LDA,
              WR3,
              WI3,
              DUMMA.asMatrix(LDU),
              LDU,
              EVECTX,
              LDU,
              N1,
              IN,
              WORK,
              IWORK,
              IWORK,
              IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DHSEIN(R)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            if (IINFO.value < 0) break;
          } else {
            // Test 11:  | HX - XW | / ( |H| |X| ulp )

            // (from inverse iteration)

            dget22('N', 'N', 'N', N, H, LDA, EVECTX, LDU, WR3, WI3, WORK,
                DUMMA(1));
            if (DUMMA[1] < ULPINV) RESULT[11] = DUMMA[1] * ANINV;
            if (DUMMA[2] > THRESH) {
              print9998(NOUNIT, 'Right', 'DHSEIN', DUMMA[2], N, JTYPE, IOLDSD);
            }
          }

          // Call DHSEIN for Left eigenvectors of H, do test 12

          NTEST = 12;
          RESULT[12] = ULPINV;
          for (var J = 1; J <= N; J++) {
            SELECT[J] = true;
          }

          dhsein('Left', 'Qr', 'Ninitv', SELECT, N, H, LDA, WR3, WI3, EVECTY,
              LDU, DUMMA.asMatrix(LDU), LDU, N1, IN, WORK, IWORK, IWORK, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DHSEIN(L)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            if (IINFO.value < 0) break;
          } else {
            // Test 12:  | YH - WY | / ( |H| |Y| ulp )

            // (from inverse iteration)

            dget22('C', 'N', 'C', N, H, LDA, EVECTY, LDU, WR3, WI3, WORK,
                DUMMA(3));
            if (DUMMA[3] < ULPINV) RESULT[12] = DUMMA[3] * ANINV;
            if (DUMMA[4] > THRESH) {
              print9998(NOUNIT, 'Left', 'DHSEIN', DUMMA[4], N, JTYPE, IOLDSD);
            }
          }

          // Call DORMHR for Right eigenvectors of A, do test 13

          NTEST = 13;
          RESULT[13] = ULPINV;

          dormhr('Left', 'No transpose', N, N, ILO, IHI, UU, LDU, TAU, EVECTX,
              LDU, WORK, NWORK, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DORMHR(R)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            if (IINFO.value < 0) break;
          } else {
            // Test 13:  | AX - XW | / ( |A| |X| ulp )

            // (from inverse iteration)

            dget22('N', 'N', 'N', N, A, LDA, EVECTX, LDU, WR3, WI3, WORK,
                DUMMA(1));
            if (DUMMA[1] < ULPINV) RESULT[13] = DUMMA[1] * ANINV;
          }

          // Call DORMHR for Left eigenvectors of A, do test 14

          NTEST = 14;
          RESULT[14] = ULPINV;

          dormhr('Left', 'No transpose', N, N, ILO, IHI, UU, LDU, TAU, EVECTY,
              LDU, WORK, NWORK, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DORMHR(L)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            if (IINFO.value < 0) break;
          } else {
            // Test 14:  | YA - WY | / ( |A| |Y| ulp )

            // (from inverse iteration)

            dget22('C', 'N', 'C', N, A, LDA, EVECTY, LDU, WR3, WI3, WORK,
                DUMMA(3));
            if (DUMMA[3] < ULPINV) RESULT[14] = DUMMA[3] * ANINV;
          }

          // Compute Left and Right Eigenvectors of A

          // Compute a Right eigenvector matrix:

          NTEST = 15;
          RESULT[15] = ULPINV;

          dlacpy(' ', N, N, UZ, LDU, EVECTR, LDU);

          dtrevc3('Right', 'Back', SELECT, N, T1, LDA, DUMMA.asMatrix(LDU), LDU,
              EVECTR, LDU, N, IN, WORK, NWORK, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DTREVC3(R,B)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            break;
          }

          // Test 15:  | AR - RW | / ( |A| |R| ulp )

          // (from Schur decomposition)

          dget22(
              'N', 'N', 'N', N, A, LDA, EVECTR, LDU, WR1, WI1, WORK, DUMMA(1));
          RESULT[15] = DUMMA[1];
          if (DUMMA[2] > THRESH) {
            print9998(NOUNIT, 'Right', 'DTREVC3', DUMMA[2], N, JTYPE, IOLDSD);
          }

          // Compute a Left eigenvector matrix:

          NTEST = 16;
          RESULT[16] = ULPINV;

          dlacpy(' ', N, N, UZ, LDU, EVECTL, LDU);

          dtrevc3('Left', 'Back', SELECT, N, T1, LDA, EVECTL, LDU,
              DUMMA.asMatrix(LDU), LDU, N, IN, WORK, NWORK, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DTREVC3(L,B)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            break;
          }

          // Test 16:  | LA - WL | / ( |A| |L| ulp )

          // (from Schur decomposition)

          dget22('Trans', 'N', 'Conj', N, A, LDA, EVECTL, LDU, WR1, WI1, WORK,
              DUMMA(3));
          RESULT[16] = DUMMA[3];
          if (DUMMA[4] > THRESH) {
            print9998(NOUNIT, 'Left', 'DTREVC3', DUMMA[4], N, JTYPE, IOLDSD);
          }

          break;
        }
        // End of Loop -- Check for RESULT[j] > THRESH

        NTESTT += NTEST;
        dlafts('DHS', N, N, JTYPE, NTEST, RESULT, IOLDSD, THRESH, NOUNIT, NERRS,
            test);
      }, skip: skip);
    }
  }

  // Summary
  dlasum('DHS', NOUNIT, NERRS.value, NTESTT);
}

void print9999(
  final Nout NOUNIT,
  final String s,
  final int info,
  final int n,
  final int ntype,
  final Array<int> iseed,
) {
  NOUNIT.println(
      ' DCHKHS: $s returned INFO=${info.i6}.\n         N=${n.i6}, JTYPE=${ntype.i6}, ISEED=(${iseed.i4(4, ',')})');
}

void print9997(
  final Nout NOUNIT,
  final String s1,
  final String s2,
  final int n,
  final int jtype,
  final Array<int> iseed,
) {
  NOUNIT.println(
      ' DCHKHS: Selected $s1 Eigenvectors from $s2 do not match other eigenvectors          N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i4(4, ',')})');
}

void print9998(
  final Nout NOUNIT,
  final String s1,
  final String s2,
  final double error,
  final int n,
  final int jtype,
  final Array<int> iseed,
) {
  NOUNIT.println(
      ' DCHKHS: $s1 Eigenvectors from $s2 incorrectly normalized.\n Bits of error=${error.g10_3},         N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i4(4, ',')})');
}
