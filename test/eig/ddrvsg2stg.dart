// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dsbgv.dart';
import 'package:dart_lapack/src/dsbgvd.dart';
import 'package:dart_lapack/src/dsbgvx.dart';
import 'package:dart_lapack/src/dspgv.dart';
import 'package:dart_lapack/src/dspgvd.dart';
import 'package:dart_lapack/src/dspgvx.dart';
import 'package:dart_lapack/src/dsygv.dart';
import 'package:dart_lapack/src/dsygv_2stage.dart';
import 'package:dart_lapack/src/dsygvd.dart';
import 'package:dart_lapack/src/dsygvx.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/range.dart';
import 'package:dart_lapack/src/xerbla.dart';

import '../matgen/dlarnd.dart';
import '../matgen/dlatmr.dart';
import '../matgen/dlatms.dart';
import '../test_driver.dart';
import 'dlafts.dart';
import 'dlasum.dart';
import 'dsgt01.dart';

void ddrvsg2stg(
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
  final int LDB,
  final Array<double> D_,
  final Array<double> D2_,
  final Matrix<double> Z_,
  final int LDZ,
  final Matrix<double> AB_,
  final Matrix<double> BB_,
  final Array<double> AP_,
  final Array<double> BP_,
  final Array<double> WORK_,
  final int NWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Array<double> RESULT_,
  final Box<int> INFO,
  final TestDriver test,
) {
  final NN = NN_.having();
  final DOTYPE = DOTYPE_.having();
  final ISEED = ISEED_.having();
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final D = D_.having();
  final D2 = D2_.having();
  final Z = Z_.having(ld: LDZ);
  final AB = AB_.having(ld: LDA);
  final BB = BB_.having(ld: LDB);
  final AP = AP_.having();
  final BP = BP_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  final RESULT = RESULT_.having();
  const ZERO = 0.0, ONE = 1.0, TEN = 10.0;
  const MAXTYP = 21;
  const KTYPE = [
    1, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8, 9, 9, 9, 9, 9, 9, //
  ];
  const KMAGN = [
    1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 1, 1, 1, 1, 1, 1, //
  ];
  const KMODE = [
    0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 4, 4, 4, 4, 4, 4, //
  ];
  // The following values are used for the half-bandwidths (ITYPE == 9)
  final KA9KB9 = [
    (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), //
    (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (1, 1), //
    (2, 1), (2, 2), (3, 1), (3, 2), (3, 3),
  ];

  // 1) Check for errors

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
  } else if (LDA <= 1 || LDA < NMAX) {
    INFO.value = -9;
  } else if (LDZ <= 1 || LDZ < NMAX) {
    INFO.value = -16;
  } else if (2 * pow(max(NMAX, 3), 2) > NWORK) {
    INFO.value = -21;
  } else if (2 * pow(max(NMAX, 3), 2) > LIWORK) {
    INFO.value = -23;
  }

  if (INFO.value != 0) {
    xerbla('DDRVSG2STG', -INFO.value);
    return;
  }

  // Quick return if possible

  if (NSIZES == 0 || NTYPES == 0) return;

  // More Important constants

  final UNFL = dlamch('Safe minimum');
  final OVFL = dlamch('Overflow');
  final ULP = dlamch('Epsilon') * dlamch('Base');
  final ULPINV = ONE / ULP;
  final RTUNFL = sqrt(UNFL);
  final RTOVFL = sqrt(OVFL);

  final ISEED2 = ISEED.copy();

  // Loop over sizes, types

  final NERRS = Box(0);

  for (final JSIZE in 1.through(NSIZES)) {
    final N = NN[JSIZE];
    final ANINV = ONE / max(1, N);
    final MTYPES = NSIZES != 1 ? min(MAXTYP, NTYPES) : min(MAXTYP + 1, NTYPES);

    for (final JTYPE in 1.through(MTYPES)) {
      final skip = !DOTYPE[JTYPE];
      test('DDRVSG2STG (SIZE = $N, TYPE = $JTYPE)', () {
        var NTEST = 0;
        final IOLDSD = ISEED.copy();
        final IINFO = Box(0), M = Box(0);
        var VL = ZERO, VU = ZERO;

        // 2)      Compute "A"
        //
        //         Control parameters:
        //
        // KMAGN  KMODE        KTYPE
        //    =1  O(1)   clustered 1  zero
        //    =2  large  clustered 2  identity
        //    =3  small  exponential  (none)
        //    =4         arithmetic   diagonal, w/ eigenvalues
        //    =5         random log   hermitian, w/ eigenvalues
        //    =6         random       (none)
        //    =7                      random diagonal
        //    =8                      random hermitian
        //    =9                      banded, w/ eigenvalues

        // Compute norm
        final ANORM = switch (KMAGN[JTYPE - 1]) {
          1 => ONE,
          2 => (RTOVFL * ULP) * ANINV,
          3 => RTUNFL * N * ULPINV,
          _ => throw UnimplementedError(),
        };

        var KA = 0, KB = 0;
        if (MTYPES <= MAXTYP) {
          final ITYPE = KTYPE[JTYPE - 1];
          final IMODE = KMODE[JTYPE - 1];
          final (KA9, KB9) = KA9KB9[JTYPE - 1];
          final COND = ULPINV;

          // Special Matrices -- Identity & Jordan block

          if (ITYPE == 1) {
            // Zero

            KA = 0;
            KB = 0;
            dlaset('Full', LDA, N, ZERO, ZERO, A, LDA);
          } else if (ITYPE == 2) {
            // Identity

            KA = 0;
            KB = 0;
            dlaset('Full', LDA, N, ZERO, ZERO, A, LDA);
            for (var JCOL = 1; JCOL <= N; JCOL++) {
              A[JCOL][JCOL] = ANORM;
            }
          } else if (ITYPE == 4) {
            // Diagonal Matrix, [Eigen]values Specified

            KA = 0;
            KB = 0;
            dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'N',
                A, LDA, WORK(N + 1), IINFO);
          } else if (ITYPE == 5) {
            // symmetric, eigenvalues specified

            KA = max(0, N - 1);
            KB = KA;
            dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, N, N, 'N',
                A, LDA, WORK(N + 1), IINFO);
          } else if (ITYPE == 7) {
            // Diagonal, random eigenvalues

            KA = 0;
            KB = 0;
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
            // symmetric, random eigenvalues

            KA = max(0, N - 1);
            KB = KA;
            final IDUMMA = Array<int>(1);
            dlatmr(
                N,
                N,
                'S',
                ISEED,
                'H',
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
            // symmetric banded, eigenvalues specified
            KA = max(0, min(N - 1, KA9));
            KB = max(0, min(N - 1, KB9));
            dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, KA, KA, 'N',
                A, LDA, WORK(N + 1), IINFO);
          } else {
            IINFO.value = 1;
          }

          test.expect(IINFO.value, 0);
          if (IINFO.value != 0) {
            _print9999(NOUNIT, 'Generator', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            return;
          }
        }

        final ABSTOL = UNFL + UNFL;
        final int IL, IU;
        if (N <= 1) {
          IL = 1;
          IU = N;
        } else {
          final TIL = 1 + ((N - 1) * dlarnd(1, ISEED2)).toInt();
          final TIU = 1 + ((N - 1) * dlarnd(1, ISEED2)).toInt();
          (IL, IU) = TIL > TIU ? (TIU, TIL) : (TIL, TIU);
        }

        // 3) Call DSYGV, DSPGV, DSBGV, SSYGVD, SSPGVD, SSBGVD,
        //    DSYGVX, DSPGVX, and DSBGVX, do tests.

        // loop over the three generalized problems
        //       IBTYPE = 1: A*x = (lambda)*B*x
        //       IBTYPE = 2: A*B*x = (lambda)*x
        //       IBTYPE = 3: B*A*x = (lambda)*x

        for (var IBTYPE = 1; IBTYPE <= 3; IBTYPE++) {
          // loop over the setting UPLO

          for (final UPLO in ['U', 'L']) {
            while (true) {
              // Generate random well-conditioned positive definite
              // matrix B, of bandwidth not greater than that of A.

              dlatms(N, N, 'U', ISEED, 'P', WORK, 5, TEN, ONE, KB, KB, UPLO, B,
                  LDB, WORK(N + 1), IINFO);

              // Test DSYGV

              NTEST++;

              dlacpy(' ', N, N, A, LDA, Z, LDZ);
              dlacpy(UPLO, N, N, B, LDB, BB, LDB);

              dsygv(
                  IBTYPE, 'V', UPLO, N, Z, LDZ, BB, LDB, D, WORK, NWORK, IINFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'DSYGV(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              dsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK,
                  RESULT(NTEST));

              // Test DSYGV_2STAGE

              NTEST++;

              dlacpy(' ', N, N, A, LDA, Z, LDZ);
              dlacpy(UPLO, N, N, B, LDB, BB, LDB);

              dsygv_2stage(IBTYPE, 'N', UPLO, N, Z, LDZ, BB, LDB, D2, WORK,
                  NWORK, IINFO);
              if (IINFO.value != 0) {
                _print9999(NOUNIT, 'DSYGV_2STAGE(V,$UPLO)', IINFO.value, N,
                    JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              // CALL DSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z,
              //              LDZ, D, WORK, RESULT[ NTEST ] )
              //
              // Do Tests | D1 - D2 | / ( |D1| ulp )
              // D1 computed using the standard 1-stage reduction as reference
              // D2 computed using the 2-stage reduction

              var TEMP1 = ZERO;
              var TEMP2 = ZERO;
              for (var J = 1; J <= N; J++) {
                TEMP1 = max(TEMP1, max(D[J].abs(), D2[J].abs()));
                TEMP2 = max(TEMP2, (D[J] - D2[J]).abs());
              }

              RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

              // Test DSYGVD

              NTEST++;

              dlacpy(' ', N, N, A, LDA, Z, LDZ);
              dlacpy(UPLO, N, N, B, LDB, BB, LDB);

              dsygvd(IBTYPE, 'V', UPLO, N, Z, LDZ, BB, LDB, D, WORK, NWORK,
                  IWORK, LIWORK, IINFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'DSYGVD(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              dsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK,
                  RESULT(NTEST));

              // Test DSYGVX

              NTEST++;

              dlacpy(' ', N, N, A, LDA, AB, LDA);
              dlacpy(UPLO, N, N, B, LDB, BB, LDB);

              dsygvx(
                  IBTYPE,
                  'V',
                  'A',
                  UPLO,
                  N,
                  AB,
                  LDA,
                  BB,
                  LDB,
                  VL,
                  VU,
                  IL,
                  IU,
                  ABSTOL,
                  M,
                  D,
                  Z,
                  LDZ,
                  WORK,
                  NWORK,
                  IWORK(N + 1),
                  IWORK,
                  IINFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'DSYGVX(V,A$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              dsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK,
                  RESULT(NTEST));

              NTEST++;

              dlacpy(' ', N, N, A, LDA, AB, LDA);
              dlacpy(UPLO, N, N, B, LDB, BB, LDB);

              // since we do not know the exact eigenvalues of this
              // eigenpair, we just set VL and VU as constants.
              // It is quite possible that there are no eigenvalues
              // in this interval.

              VL = ZERO;
              VU = ANORM;
              dsygvx(
                  IBTYPE,
                  'V',
                  'V',
                  UPLO,
                  N,
                  AB,
                  LDA,
                  BB,
                  LDB,
                  VL,
                  VU,
                  IL,
                  IU,
                  ABSTOL,
                  M,
                  D,
                  Z,
                  LDZ,
                  WORK,
                  NWORK,
                  IWORK(N + 1),
                  IWORK,
                  IINFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'DSYGVX(V,V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              dsgt01(IBTYPE, UPLO, N, M.value, A, LDA, B, LDB, Z, LDZ, D, WORK,
                  RESULT(NTEST));

              NTEST++;

              dlacpy(' ', N, N, A, LDA, AB, LDA);
              dlacpy(UPLO, N, N, B, LDB, BB, LDB);

              dsygvx(
                  IBTYPE,
                  'V',
                  'I',
                  UPLO,
                  N,
                  AB,
                  LDA,
                  BB,
                  LDB,
                  VL,
                  VU,
                  IL,
                  IU,
                  ABSTOL,
                  M,
                  D,
                  Z,
                  LDZ,
                  WORK,
                  NWORK,
                  IWORK(N + 1),
                  IWORK,
                  IINFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'DSYGVX(V,I,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              dsgt01(IBTYPE, UPLO, N, M.value, A, LDA, B, LDB, Z, LDZ, D, WORK,
                  RESULT(NTEST));

              break;
            }

            // Test DSPGV
            while (true) {
              NTEST++;

              // Copy the matrices into packed storage.

              if (lsame(UPLO, 'U')) {
                for (var J = 1, IJ = 1; J <= N; J++) {
                  for (var I = 1; I <= J; I++) {
                    AP[IJ] = A[I][J];
                    BP[IJ] = B[I][J];
                    IJ++;
                  }
                }
              } else {
                for (var J = 1, IJ = 1; J <= N; J++) {
                  for (var I = J; I <= N; I++) {
                    AP[IJ] = A[I][J];
                    BP[IJ] = B[I][J];
                    IJ++;
                  }
                }
              }

              dspgv(IBTYPE, 'V', UPLO, N, AP, BP, D, Z, LDZ, WORK, IINFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'DSPGV(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              dsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK,
                  RESULT(NTEST));

              // Test DSPGVD

              NTEST++;

              // Copy the matrices into packed storage.

              if (lsame(UPLO, 'U')) {
                for (var J = 1, IJ = 1; J <= N; J++) {
                  for (var I = 1; I <= J; I++) {
                    AP[IJ] = A[I][J];
                    BP[IJ] = B[I][J];
                    IJ++;
                  }
                }
              } else {
                for (var J = 1, IJ = 1; J <= N; J++) {
                  for (var I = J; I <= N; I++) {
                    AP[IJ] = A[I][J];
                    BP[IJ] = B[I][J];
                    IJ++;
                  }
                }
              }

              dspgvd(IBTYPE, 'V', UPLO, N, AP, BP, D, Z, LDZ, WORK, NWORK,
                  IWORK, LIWORK, IINFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'DSPGVD(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              dsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK,
                  RESULT(NTEST));

              // Test DSPGVX

              NTEST++;

              // Copy the matrices into packed storage.

              if (lsame(UPLO, 'U')) {
                for (var J = 1, IJ = 1; J <= N; J++) {
                  for (var I = 1; I <= J; I++) {
                    AP[IJ] = A[I][J];
                    BP[IJ] = B[I][J];
                    IJ++;
                  }
                }
              } else {
                for (var J = 1, IJ = 1; J <= N; J++) {
                  for (var I = J; I <= N; I++) {
                    AP[IJ] = A[I][J];
                    BP[IJ] = B[I][J];
                    IJ++;
                  }
                }
              }

              dspgvx(IBTYPE, 'V', 'A', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL,
                  M, D, Z, LDZ, WORK, IWORK(N + 1), IWORK, INFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'DSPGVX(V,A$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              dsgt01(IBTYPE, UPLO, N, M.value, A, LDA, B, LDB, Z, LDZ, D, WORK,
                  RESULT(NTEST));

              NTEST++;

              // Copy the matrices into packed storage.

              if (lsame(UPLO, 'U')) {
                for (var J = 1, IJ = 1; J <= N; J++) {
                  for (var I = 1; I <= J; I++) {
                    AP[IJ] = A[I][J];
                    BP[IJ] = B[I][J];
                    IJ++;
                  }
                }
              } else {
                for (var J = 1, IJ = 1; J <= N; J++) {
                  for (var I = J; I <= N; I++) {
                    AP[IJ] = A[I][J];
                    BP[IJ] = B[I][J];
                    IJ++;
                  }
                }
              }

              VL = ZERO;
              VU = ANORM;
              dspgvx(IBTYPE, 'V', 'V', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL,
                  M, D, Z, LDZ, WORK, IWORK(N + 1), IWORK, INFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'DSPGVX(V,V$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              dsgt01(IBTYPE, UPLO, N, M.value, A, LDA, B, LDB, Z, LDZ, D, WORK,
                  RESULT(NTEST));

              NTEST++;

              // Copy the matrices into packed storage.

              if (lsame(UPLO, 'U')) {
                for (var J = 1, IJ = 1; J <= N; J++) {
                  for (var I = 1; I <= J; I++) {
                    AP[IJ] = A[I][J];
                    BP[IJ] = B[I][J];
                    IJ++;
                  }
                }
              } else {
                for (var J = 1, IJ = 1; J <= N; J++) {
                  for (var I = J; I <= N; I++) {
                    AP[IJ] = A[I][J];
                    BP[IJ] = B[I][J];
                    IJ++;
                  }
                }
              }

              dspgvx(IBTYPE, 'V', 'I', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL,
                  M, D, Z, LDZ, WORK, IWORK(N + 1), IWORK, INFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'DSPGVX(V,I$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              dsgt01(IBTYPE, UPLO, N, M.value, A, LDA, B, LDB, Z, LDZ, D, WORK,
                  RESULT(NTEST));

              break;
            }

            while (IBTYPE == 1) {
              // TEST DSBGV

              NTEST++;

              // Copy the matrices into band storage.

              if (lsame(UPLO, 'U')) {
                for (var J = 1; J <= N; J++) {
                  for (var I = max(1, J - KA); I <= J; I++) {
                    AB[KA + 1 + I - J][J] = A[I][J];
                  }
                  for (var I = max(1, J - KB); I <= J; I++) {
                    BB[KB + 1 + I - J][J] = B[I][J];
                  }
                }
              } else {
                for (var J = 1; J <= N; J++) {
                  for (var I = J; I <= min(N, J + KA); I++) {
                    AB[1 + I - J][J] = A[I][J];
                  }
                  for (var I = J; I <= min(N, J + KB); I++) {
                    BB[1 + I - J][J] = B[I][J];
                  }
                }
              }

              dsbgv('V', UPLO, N, KA, KB, AB, LDA, BB, LDB, D, Z, LDZ, WORK,
                  IINFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'DSBGV(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              dsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK,
                  RESULT(NTEST));

              // TEST DSBGVD

              NTEST++;

              // Copy the matrices into band storage.

              if (lsame(UPLO, 'U')) {
                for (var J = 1; J <= N; J++) {
                  for (var I = max(1, J - KA); I <= J; I++) {
                    AB[KA + 1 + I - J][J] = A[I][J];
                  }
                  for (var I = max(1, J - KB); I <= J; I++) {
                    BB[KB + 1 + I - J][J] = B[I][J];
                  }
                }
              } else {
                for (var J = 1; J <= N; J++) {
                  for (var I = J; I <= min(N, J + KA); I++) {
                    AB[1 + I - J][J] = A[I][J];
                  }
                  for (var I = J; I <= min(N, J + KB); I++) {
                    BB[1 + I - J][J] = B[I][J];
                  }
                }
              }

              dsbgvd('V', UPLO, N, KA, KB, AB, LDA, BB, LDB, D, Z, LDZ, WORK,
                  NWORK, IWORK, LIWORK, IINFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'DSBGVD(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              dsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK,
                  RESULT(NTEST));

              // Test DSBGVX

              NTEST++;

              // Copy the matrices into band storage.

              if (lsame(UPLO, 'U')) {
                for (var J = 1; J <= N; J++) {
                  for (var I = max(1, J - KA); I <= J; I++) {
                    AB[KA + 1 + I - J][J] = A[I][J];
                  }
                  for (var I = max(1, J - KB); I <= J; I++) {
                    BB[KB + 1 + I - J][J] = B[I][J];
                  }
                }
              } else {
                for (var J = 1; J <= N; J++) {
                  for (var I = J; I <= min(N, J + KA); I++) {
                    AB[1 + I - J][J] = A[I][J];
                  }
                  for (var I = J; I <= min(N, J + KB); I++) {
                    BB[1 + I - J][J] = B[I][J];
                  }
                }
              }

              dsbgvx(
                  'V',
                  'A',
                  UPLO,
                  N,
                  KA,
                  KB,
                  AB,
                  LDA,
                  BB,
                  LDB,
                  BP.asMatrix(max(1, N)),
                  max(1, N),
                  VL,
                  VU,
                  IL,
                  IU,
                  ABSTOL,
                  M,
                  D,
                  Z,
                  LDZ,
                  WORK,
                  IWORK(N + 1),
                  IWORK,
                  IINFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'DSBGVX(V,A$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              dsgt01(IBTYPE, UPLO, N, M.value, A, LDA, B, LDB, Z, LDZ, D, WORK,
                  RESULT(NTEST));

              NTEST++;

              // Copy the matrices into band storage.

              if (lsame(UPLO, 'U')) {
                for (var J = 1; J <= N; J++) {
                  for (var I = max(1, J - KA); I <= J; I++) {
                    AB[KA + 1 + I - J][J] = A[I][J];
                  }
                  for (var I = max(1, J - KB); I <= J; I++) {
                    BB[KB + 1 + I - J][J] = B[I][J];
                  }
                }
              } else {
                for (var J = 1; J <= N; J++) {
                  for (var I = J; I <= min(N, J + KA); I++) {
                    AB[1 + I - J][J] = A[I][J];
                  }
                  for (var I = J; I <= min(N, J + KB); I++) {
                    BB[1 + I - J][J] = B[I][J];
                  }
                }
              }

              VL = ZERO;
              VU = ANORM;
              dsbgvx(
                  'V',
                  'V',
                  UPLO,
                  N,
                  KA,
                  KB,
                  AB,
                  LDA,
                  BB,
                  LDB,
                  BP.asMatrix(max(1, N)),
                  max(1, N),
                  VL,
                  VU,
                  IL,
                  IU,
                  ABSTOL,
                  M,
                  D,
                  Z,
                  LDZ,
                  WORK,
                  IWORK(N + 1),
                  IWORK,
                  IINFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'DSBGVX(V,V$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              dsgt01(IBTYPE, UPLO, N, M.value, A, LDA, B, LDB, Z, LDZ, D, WORK,
                  RESULT(NTEST));

              NTEST++;

              // Copy the matrices into band storage.

              if (lsame(UPLO, 'U')) {
                for (var J = 1; J <= N; J++) {
                  for (var I = max(1, J - KA); I <= J; I++) {
                    AB[KA + 1 + I - J][J] = A[I][J];
                  }
                  for (var I = max(1, J - KB); I <= J; I++) {
                    BB[KB + 1 + I - J][J] = B[I][J];
                  }
                }
              } else {
                for (var J = 1; J <= N; J++) {
                  for (var I = J; I <= min(N, J + KA); I++) {
                    AB[1 + I - J][J] = A[I][J];
                  }
                  for (var I = J; I <= min(N, J + KB); I++) {
                    BB[1 + I - J][J] = B[I][J];
                  }
                }
              }

              dsbgvx(
                  'V',
                  'I',
                  UPLO,
                  N,
                  KA,
                  KB,
                  AB,
                  LDA,
                  BB,
                  LDB,
                  BP.asMatrix(max(1, N)),
                  max(1, N),
                  VL,
                  VU,
                  IL,
                  IU,
                  ABSTOL,
                  M,
                  D,
                  Z,
                  LDZ,
                  WORK,
                  IWORK(N + 1),
                  IWORK,
                  IINFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'DSBGVX(V,I$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              dsgt01(IBTYPE, UPLO, N, M.value, A, LDA, B, LDB, Z, LDZ, D, WORK,
                  RESULT(NTEST));

              break;
            }
          }
        }

        // End of Loop -- Check for RESULT[j] > THRESH

        NTESTT += NTEST;
        dlafts('DSG', N, N, JTYPE, NTEST, RESULT, IOLDSD, THRESH, NOUNIT, NERRS,
            test);
      }, skip: skip);
    }
  }

  // Summary
  dlasum('DSG', NOUNIT, NERRS.value, NTESTT);
}

void _print9999(
  final Nout nout,
  final String s,
  final int info,
  final int n,
  final int jtype,
  final Array<int> iseed,
) {
  nout.println(
      ' DDRVSG2STG: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}
