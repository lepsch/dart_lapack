// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhbgv.dart';
import 'package:lapack/src/zhbgvd.dart';
import 'package:lapack/src/zhbgvx.dart';
import 'package:lapack/src/zhegv.dart';
import 'package:lapack/src/zhegvd.dart';
import 'package:lapack/src/zhegvx.dart';
import 'package:lapack/src/zhpgv.dart';
import 'package:lapack/src/zhpgvd.dart';
import 'package:lapack/src/zhpgvx.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlaset.dart';

import '../matgen/dlarnd.dart';
import '../matgen/zlatmr.dart';
import '../matgen/zlatms.dart';
import 'dlafts.dart';
import 'dlasum.dart';
import 'zsgt01.dart';

void zdrvsg(
  final int NSIZES,
  final Array<int> NN_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final Array<int> ISEED_,
  final double THRESH,
  final Nout NOUNIT,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<double> D_,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Matrix<Complex> AB_,
  final Matrix<Complex> BB_,
  final Array<Complex> AP_,
  final Array<Complex> BP_,
  final Array<Complex> WORK_,
  final int NWORK,
  final Array<double> RWORK_,
  final int LRWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Array<double> RESULT_,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final NN = NN_.having();
  final DOTYPE = DOTYPE_.having();
  final ISEED = ISEED_.having(length: 4);
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final Z = Z_.having(ld: LDZ);
  final AB = AB_.having(ld: LDA);
  final BB = BB_.having(ld: LDB);
  final AP = AP_.having();
  final BP = BP_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  final D = D_.having();
  final RESULT = RESULT_.having();
  const ZERO = 0.0, ONE = 1.0, TEN = 10.0;
  const MAXTYP = 21;
  bool BADNN;
  String UPLO = '';
  int I,
      IBTYPE,
      IBUPLO,
      IJ,
      IL,
      IMODE,
      ITEMP,
      ITYPE,
      IU,
      J,
      JCOL,
      JSIZE,
      JTYPE,
      KA = 0,
      KA9,
      KB = 0,
      KB9,
      MTYPES,
      N,
      NMAX,
      NTEST,
      NTESTT;
  double ABSTOL,
      ANINV,
      ANORM = 0,
      COND,
      OVFL,
      RTOVFL,
      RTUNFL,
      ULP,
      ULPINV,
      UNFL,
      VL = 0,
      VU = 0;
  final IDUMMA = Array<int>(1), IOLDSD = Array<int>(4), ISEED2 = Array<int>(4);
  final KTYPE = Array.fromList(
      [1, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8, 9, 9, 9, 9, 9, 9]);
  final KMAGN = Array.fromList(
      [1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 1, 1, 1, 1, 1, 1]);
  final KMODE = Array.fromList(
      [0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 4, 4, 4, 4, 4, 4]);
  final IINFO = Box(0), M = Box(0), NERRS = Box(0);

  // 1)      Check for errors

  NTESTT = 0;
  INFO.value = 0;

  BADNN = false;
  NMAX = 0;
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
  } else if (LDA <= 1 || LDA < NMAX) {
    INFO.value = -9;
  } else if (LDZ <= 1 || LDZ < NMAX) {
    INFO.value = -16;
  } else if (2 * pow(max(NMAX, 2), 2) > NWORK) {
    INFO.value = -21;
  } else if (2 * pow(max(NMAX, 2), 2) > LRWORK) {
    INFO.value = -23;
  } else if (2 * pow(max(NMAX, 2), 2) > LIWORK) {
    INFO.value = -25;
  }

  if (INFO.value != 0) {
    xerbla('ZDRVSG', -INFO.value);
    return;
  }

  // Quick return if possible

  if (NSIZES == 0 || NTYPES == 0) return;

  // More Important constants

  UNFL = dlamch('Safe minimum');
  OVFL = dlamch('Overflow');
  ULP = dlamch('Epsilon') * dlamch('Base');
  ULPINV = ONE / ULP;
  RTUNFL = sqrt(UNFL);
  RTOVFL = sqrt(OVFL);

  for (I = 1; I <= 4; I++) {
    ISEED2[I] = ISEED[I];
  }

  // Loop over sizes, types

  NERRS.value = 0;

  for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) {
    N = NN[JSIZE];
    ANINV = ONE / max(1, N);

    if (NSIZES != 1) {
      MTYPES = min(MAXTYP, NTYPES);
    } else {
      MTYPES = min(MAXTYP + 1, NTYPES);
    }

    KA9 = 0;
    KB9 = 0;
    for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
      if (!DOTYPE[JTYPE]) continue;
      NTEST = 0;

      for (J = 1; J <= 4; J++) {
        IOLDSD[J] = ISEED[J];
      }

      // 2)      Compute "A"
      //
      //         Control parameters:
      //
      //     KMAGN  KMODE        KTYPE
      // =1  O(1)   clustered 1  zero
      // =2  large  clustered 2  identity
      // =3  small  exponential  (none)
      // =4         arithmetic   diagonal, w/ eigenvalues
      // =5         random log   hermitian, w/ eigenvalues
      // =6         random       (none)
      // =7                      random diagonal
      // =8                      random hermitian
      // =9                      banded, w/ eigenvalues

      if (MTYPES <= MAXTYP) {
        ITYPE = KTYPE[JTYPE];
        IMODE = KMODE[JTYPE];

        // Compute norm

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

        IINFO.value = 0;
        COND = ULPINV;

        // Special Matrices -- Identity & Jordan block

        if (ITYPE == 1) {
          // Zero

          KA = 0;
          KB = 0;
          zlaset('Full', LDA, N, Complex.zero, Complex.zero, A, LDA);
        } else if (ITYPE == 2) {
          // Identity

          KA = 0;
          KB = 0;
          zlaset('Full', LDA, N, Complex.zero, Complex.zero, A, LDA);
          for (JCOL = 1; JCOL <= N; JCOL++) {
            A[JCOL][JCOL] = ANORM.toComplex();
          }
        } else if (ITYPE == 4) {
          // Diagonal Matrix, [Eigen]values Specified

          KA = 0;
          KB = 0;
          zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, 0, 0, 'N', A,
              LDA, WORK, IINFO);
        } else if (ITYPE == 5) {
          // Hermitian, eigenvalues specified

          KA = max(0, N - 1);
          KB = KA;
          zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, N, N, 'N', A,
              LDA, WORK, IINFO);
        } else if (ITYPE == 7) {
          // Diagonal, random eigenvalues

          KA = 0;
          KB = 0;
          zlatmr(
              N,
              N,
              'S',
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
          // Hermitian, random eigenvalues

          KA = max(0, N - 1);
          KB = KA;
          zlatmr(
              N,
              N,
              'S',
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
              IWORK,
              IINFO);
        } else if (ITYPE == 9) {
          // Hermitian banded, eigenvalues specified

          // The following values are used for the half-bandwidths:
          //
          //   ka = 1   kb = 1
          //   ka = 2   kb = 1
          //   ka = 2   kb = 2
          //   ka = 3   kb = 1
          //   ka = 3   kb = 2
          //   ka = 3   kb = 3

          KB9++;
          if (KB9 > KA9) {
            KA9++;
            KB9 = 1;
          }
          KA = max(0, min(N - 1, KA9));
          KB = max(0, min(N - 1, KB9));
          zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, KA, KA, 'N',
              A, LDA, WORK, IINFO);
        } else {
          IINFO.value = 1;
        }

        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'Generator', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }
      }

      ABSTOL = UNFL + UNFL;
      if (N <= 1) {
        IL = 1;
        IU = N;
      } else {
        IL = 1 + ((N - 1) * dlarnd(1, ISEED2)).toInt();
        IU = 1 + ((N - 1) * dlarnd(1, ISEED2)).toInt();
        if (IL > IU) {
          ITEMP = IL;
          IL = IU;
          IU = ITEMP;
        }
      }

      // 3) Call ZHEGV, ZHPGV, ZHBGV, CHEGVD, CHPGVD, CHBGVD,
      //    ZHEGVX, ZHPGVX and ZHBGVX, do tests.

      // loop over the three generalized problems
      //       IBTYPE = 1: A*x = (lambda)*B*x
      //       IBTYPE = 2: A*B*x = (lambda)*x
      //       IBTYPE = 3: B*A*x = (lambda)*x

      for (IBTYPE = 1; IBTYPE <= 3; IBTYPE++) {
        // loop over the setting UPLO

        for (IBUPLO = 1; IBUPLO <= 2; IBUPLO++) {
          if (IBUPLO == 1) UPLO = 'U';
          if (IBUPLO == 2) UPLO = 'L';

          // Generate random well-conditioned positive definite
          // matrix B, of bandwidth not greater than that of A.

          zlatms(N, N, 'U', ISEED, 'P', RWORK, 5, TEN, ONE, KB, KB, UPLO, B,
              LDB, WORK(N + 1), IINFO);

          // Test ZHEGV

          NTEST++;

          zlacpy(' ', N, N, A, LDA, Z, LDZ);
          zlacpy(UPLO, N, N, B, LDB, BB, LDB);

          while (true) {
            zhegv(IBTYPE, 'V', UPLO, N, Z, LDZ, BB, LDB, D, WORK, NWORK, RWORK,
                IINFO);
            if (IINFO.value != 0) {
              _print9999(
                  NOUNIT, 'ZHEGV(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[NTEST] = ULPINV;
                break;
              }
            }

            // Do Test

            zsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK,
                RESULT(NTEST));

            // Test ZHEGVD

            NTEST++;

            zlacpy(' ', N, N, A, LDA, Z, LDZ);
            zlacpy(UPLO, N, N, B, LDB, BB, LDB);

            zhegvd(IBTYPE, 'V', UPLO, N, Z, LDZ, BB, LDB, D, WORK, NWORK, RWORK,
                LRWORK, IWORK, LIWORK, IINFO);
            if (IINFO.value != 0) {
              _print9999(
                  NOUNIT, 'ZHEGVD(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[NTEST] = ULPINV;
                break;
              }
            }

            // Do Test

            zsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK,
                RESULT(NTEST));

            // Test ZHEGVX

            NTEST++;

            zlacpy(' ', N, N, A, LDA, AB, LDA);
            zlacpy(UPLO, N, N, B, LDB, BB, LDB);

            zhegvx(
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
                RWORK,
                IWORK(N + 1),
                IWORK,
                IINFO);
            if (IINFO.value != 0) {
              _print9999(
                  NOUNIT, 'ZHEGVX(V,A$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[NTEST] = ULPINV;
                break;
              }
            }

            // Do Test

            zsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK,
                RESULT(NTEST));

            NTEST++;

            zlacpy(' ', N, N, A, LDA, AB, LDA);
            zlacpy(UPLO, N, N, B, LDB, BB, LDB);

            // since we do not know the exact eigenvalues of this
            // eigenpair, we just set VL and VU as constants.
            // It is quite possible that there are no eigenvalues
            // in this interval.

            VL = ZERO;
            VU = ANORM;
            zhegvx(
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
                RWORK,
                IWORK(N + 1),
                IWORK,
                IINFO);
            if (IINFO.value != 0) {
              _print9999(
                  NOUNIT, 'ZHEGVX(V,V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[NTEST] = ULPINV;
                break;
              }
            }

            // Do Test

            zsgt01(IBTYPE, UPLO, N, M.value, A, LDA, B, LDB, Z, LDZ, D, WORK,
                RWORK, RESULT(NTEST));

            NTEST++;

            zlacpy(' ', N, N, A, LDA, AB, LDA);
            zlacpy(UPLO, N, N, B, LDB, BB, LDB);

            zhegvx(
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
                RWORK,
                IWORK(N + 1),
                IWORK,
                IINFO);
            if (IINFO.value != 0) {
              _print9999(
                  NOUNIT, 'ZHEGVX(V,I,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[NTEST] = ULPINV;
                break;
              }
            }

            // Do Test

            zsgt01(IBTYPE, UPLO, N, M.value, A, LDA, B, LDB, Z, LDZ, D, WORK,
                RWORK, RESULT(NTEST));

            break;
          }

          // Test ZHPGV

          NTEST++;

          // Copy the matrices into packed storage.

          if (lsame(UPLO, 'U')) {
            IJ = 1;
            for (J = 1; J <= N; J++) {
              for (I = 1; I <= J; I++) {
                AP[IJ] = A[I][J];
                BP[IJ] = B[I][J];
                IJ++;
              }
            }
          } else {
            IJ = 1;
            for (J = 1; J <= N; J++) {
              for (I = J; I <= N; I++) {
                AP[IJ] = A[I][J];
                BP[IJ] = B[I][J];
                IJ++;
              }
            }
          }

          while (true) {
            zhpgv(IBTYPE, 'V', UPLO, N, AP, BP, D, Z, LDZ, WORK, RWORK, IINFO);
            if (IINFO.value != 0) {
              _print9999(
                  NOUNIT, 'ZHPGV(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[NTEST] = ULPINV;
                break;
              }
            }

            // Do Test

            zsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK,
                RESULT(NTEST));

            // Test ZHPGVD

            NTEST++;

            // Copy the matrices into packed storage.

            if (lsame(UPLO, 'U')) {
              IJ = 1;
              for (J = 1; J <= N; J++) {
                for (I = 1; I <= J; I++) {
                  AP[IJ] = A[I][J];
                  BP[IJ] = B[I][J];
                  IJ++;
                }
              }
            } else {
              IJ = 1;
              for (J = 1; J <= N; J++) {
                for (I = J; I <= N; I++) {
                  AP[IJ] = A[I][J];
                  BP[IJ] = B[I][J];
                  IJ++;
                }
              }
            }

            zhpgvd(IBTYPE, 'V', UPLO, N, AP, BP, D, Z, LDZ, WORK, NWORK, RWORK,
                LRWORK, IWORK, LIWORK, IINFO);
            if (IINFO.value != 0) {
              _print9999(
                  NOUNIT, 'ZHPGVD(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[NTEST] = ULPINV;
                break;
              }
            }

            // Do Test

            zsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK,
                RESULT(NTEST));

            // Test ZHPGVX

            NTEST++;

            // Copy the matrices into packed storage.

            if (lsame(UPLO, 'U')) {
              IJ = 1;
              for (J = 1; J <= N; J++) {
                for (I = 1; I <= J; I++) {
                  AP[IJ] = A[I][J];
                  BP[IJ] = B[I][J];
                  IJ++;
                }
              }
            } else {
              IJ = 1;
              for (J = 1; J <= N; J++) {
                for (I = J; I <= N; I++) {
                  AP[IJ] = A[I][J];
                  BP[IJ] = B[I][J];
                  IJ++;
                }
              }
            }

            zhpgvx(IBTYPE, 'V', 'A', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M,
                D, Z, LDZ, WORK, RWORK, IWORK(N + 1), IWORK, IINFO);
            if (IINFO.value != 0) {
              _print9999(
                  NOUNIT, 'ZHPGVX(V,A$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[NTEST] = ULPINV;
                break;
              }
            }

            // Do Test

            zsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK,
                RESULT(NTEST));

            NTEST++;

            // Copy the matrices into packed storage.

            if (lsame(UPLO, 'U')) {
              IJ = 1;
              for (J = 1; J <= N; J++) {
                for (I = 1; I <= J; I++) {
                  AP[IJ] = A[I][J];
                  BP[IJ] = B[I][J];
                  IJ++;
                }
              }
            } else {
              IJ = 1;
              for (J = 1; J <= N; J++) {
                for (I = J; I <= N; I++) {
                  AP[IJ] = A[I][J];
                  BP[IJ] = B[I][J];
                  IJ++;
                }
              }
            }

            VL = ZERO;
            VU = ANORM;
            zhpgvx(IBTYPE, 'V', 'V', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M,
                D, Z, LDZ, WORK, RWORK, IWORK(N + 1), IWORK, IINFO);
            if (IINFO.value != 0) {
              _print9999(
                  NOUNIT, 'ZHPGVX(V,V$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[NTEST] = ULPINV;
                break;
              }
            }

            // Do Test

            zsgt01(IBTYPE, UPLO, N, M.value, A, LDA, B, LDB, Z, LDZ, D, WORK,
                RWORK, RESULT(NTEST));

            NTEST++;

            // Copy the matrices into packed storage.

            if (lsame(UPLO, 'U')) {
              IJ = 1;
              for (J = 1; J <= N; J++) {
                for (I = 1; I <= J; I++) {
                  AP[IJ] = A[I][J];
                  BP[IJ] = B[I][J];
                  IJ++;
                }
              }
            } else {
              IJ = 1;
              for (J = 1; J <= N; J++) {
                for (I = J; I <= N; I++) {
                  AP[IJ] = A[I][J];
                  BP[IJ] = B[I][J];
                  IJ++;
                }
              }
            }

            zhpgvx(IBTYPE, 'V', 'I', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M,
                D, Z, LDZ, WORK, RWORK, IWORK(N + 1), IWORK, IINFO);
            if (IINFO.value != 0) {
              _print9999(
                  NOUNIT, 'ZHPGVX(V,I$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[NTEST] = ULPINV;
                break;
              }
            }

            // Do Test

            zsgt01(IBTYPE, UPLO, N, M.value, A, LDA, B, LDB, Z, LDZ, D, WORK,
                RWORK, RESULT(NTEST));

            break;
          }

          if (IBTYPE == 1) {
            // TEST ZHBGV

            NTEST++;

            // Copy the matrices into band storage.

            if (lsame(UPLO, 'U')) {
              for (J = 1; J <= N; J++) {
                for (I = max(1, J - KA); I <= J; I++) {
                  AB[KA + 1 + I - J][J] = A[I][J];
                }
                for (I = max(1, J - KB); I <= J; I++) {
                  BB[KB + 1 + I - J][J] = B[I][J];
                }
              }
            } else {
              for (J = 1; J <= N; J++) {
                for (I = J; I <= min(N, J + KA); I++) {
                  AB[1 + I - J][J] = A[I][J];
                }
                for (I = J; I <= min(N, J + KB); I++) {
                  BB[1 + I - J][J] = B[I][J];
                }
              }
            }

            while (true) {
              zhbgv('V', UPLO, N, KA, KB, AB, LDA, BB, LDB, D, Z, LDZ, WORK,
                  RWORK, IINFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'ZHBGV(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              zsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK,
                  RESULT(NTEST));

              // TEST ZHBGVD

              NTEST++;

              // Copy the matrices into band storage.

              if (lsame(UPLO, 'U')) {
                for (J = 1; J <= N; J++) {
                  for (I = max(1, J - KA); I <= J; I++) {
                    AB[KA + 1 + I - J][J] = A[I][J];
                  }
                  for (I = max(1, J - KB); I <= J; I++) {
                    BB[KB + 1 + I - J][J] = B[I][J];
                  }
                }
              } else {
                for (J = 1; J <= N; J++) {
                  for (I = J; I <= min(N, J + KA); I++) {
                    AB[1 + I - J][J] = A[I][J];
                  }
                  for (I = J; I <= min(N, J + KB); I++) {
                    BB[1 + I - J][J] = B[I][J];
                  }
                }
              }

              zhbgvd('V', UPLO, N, KA, KB, AB, LDA, BB, LDB, D, Z, LDZ, WORK,
                  NWORK, RWORK, LRWORK, IWORK, LIWORK, IINFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'ZHBGVD(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              zsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK,
                  RESULT(NTEST));

              // Test ZHBGVX

              NTEST++;

              // Copy the matrices into band storage.

              if (lsame(UPLO, 'U')) {
                for (J = 1; J <= N; J++) {
                  for (I = max(1, J - KA); I <= J; I++) {
                    AB[KA + 1 + I - J][J] = A[I][J];
                  }
                  for (I = max(1, J - KB); I <= J; I++) {
                    BB[KB + 1 + I - J][J] = B[I][J];
                  }
                }
              } else {
                for (J = 1; J <= N; J++) {
                  for (I = J; I <= min(N, J + KA); I++) {
                    AB[1 + I - J][J] = A[I][J];
                  }
                  for (I = J; I <= min(N, J + KB); I++) {
                    BB[1 + I - J][J] = B[I][J];
                  }
                }
              }

              zhbgvx(
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
                  BP.asMatrix(),
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
                  RWORK,
                  IWORK(N + 1),
                  IWORK,
                  IINFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'ZHBGVX(V,A$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              zsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK,
                  RESULT(NTEST));

              NTEST++;

              // Copy the matrices into band storage.

              if (lsame(UPLO, 'U')) {
                for (J = 1; J <= N; J++) {
                  for (I = max(1, J - KA); I <= J; I++) {
                    AB[KA + 1 + I - J][J] = A[I][J];
                  }
                  for (I = max(1, J - KB); I <= J; I++) {
                    BB[KB + 1 + I - J][J] = B[I][J];
                  }
                }
              } else {
                for (J = 1; J <= N; J++) {
                  for (I = J; I <= min(N, J + KA); I++) {
                    AB[1 + I - J][J] = A[I][J];
                  }
                  for (I = J; I <= min(N, J + KB); I++) {
                    BB[1 + I - J][J] = B[I][J];
                  }
                }
              }

              VL = ZERO;
              VU = ANORM;
              zhbgvx(
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
                  BP.asMatrix(),
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
                  RWORK,
                  IWORK(N + 1),
                  IWORK,
                  IINFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'ZHBGVX(V,V$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              zsgt01(IBTYPE, UPLO, N, M.value, A, LDA, B, LDB, Z, LDZ, D, WORK,
                  RWORK, RESULT(NTEST));

              NTEST++;

              // Copy the matrices into band storage.

              if (lsame(UPLO, 'U')) {
                for (J = 1; J <= N; J++) {
                  for (I = max(1, J - KA); I <= J; I++) {
                    AB[KA + 1 + I - J][J] = A[I][J];
                  }
                  for (I = max(1, J - KB); I <= J; I++) {
                    BB[KB + 1 + I - J][J] = B[I][J];
                  }
                }
              } else {
                for (J = 1; J <= N; J++) {
                  for (I = J; I <= min(N, J + KA); I++) {
                    AB[1 + I - J][J] = A[I][J];
                  }
                  for (I = J; I <= min(N, J + KB); I++) {
                    BB[1 + I - J][J] = B[I][J];
                  }
                }
              }

              zhbgvx(
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
                  BP.asMatrix(),
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
                  RWORK,
                  IWORK(N + 1),
                  IWORK,
                  IINFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'ZHBGVX(V,I$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              zsgt01(IBTYPE, UPLO, N, M.value, A, LDA, B, LDB, Z, LDZ, D, WORK,
                  RWORK, RESULT(NTEST));
            }
          }
        }
      }

      // End of Loop -- Check for RESULT(j) > THRESH

      NTESTT += NTEST;
      dlafts('ZSG', N, N, JTYPE, NTEST, RESULT, IOLDSD, THRESH, NOUNIT, NERRS);
    }
  }

  // Summary

  dlasum('ZSG', NOUNIT, NERRS.value, NTESTT);
}

void _print9999(
  Nout nout,
  String s,
  int info,
  int n,
  int jtype,
  Array<int> iseed,
) {
  nout.println(
      ' ZDRVSG: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}
