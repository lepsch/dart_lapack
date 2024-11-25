// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlanhe.dart';
import 'package:dart_lapack/src/zpftrf.dart';
import 'package:dart_lapack/src/zpftri.dart';
import 'package:dart_lapack/src/zpftrs.dart';
import 'package:dart_lapack/src/zpotrf.dart';
import 'package:dart_lapack/src/zpotri.dart';
import 'package:dart_lapack/src/ztfttr.dart';
import 'package:dart_lapack/src/ztrttf.dart';

import '../matgen/zlatms.dart';
import 'aladhd.dart';
import 'alaerh.dart';
import 'alasvm.dart';
import 'common.dart';
import 'zget04.dart';
import 'zlaipd.dart';
import 'zlarhs.dart';
import 'zlatb4.dart';
import 'zpot01.dart';
import 'zpot02.dart';
import 'zpot03.dart';

void zdrvrfp(
  final Nout NOUT,
  final int NN,
  final Array<int> NVAL_,
  final int NNS,
  final Array<int> NSVAL_,
  final int NNT,
  final Array<int> NTVAL,
  final double THRESH,
  final Array<Complex> A_,
  final Array<Complex> ASAV_,
  final Array<Complex> AFAC_,
  final Array<Complex> AINV_,
  final Array<Complex> B_,
  final Array<Complex> BSAV_,
  final Array<Complex> XACT_,
  final Array<Complex> X_,
  final Array<Complex> ARF_,
  final Array<Complex> ARFINV_,
  final Array<Complex> Z_WORK_ZLATMS_,
  final Array<Complex> Z_WORK_ZPOT02_,
  final Array<Complex> Z_WORK_ZPOT03_,
  final Array<double> D_WORK_ZLATMS_,
  final Array<double> D_WORK_ZLANHE_,
  final Array<double> D_WORK_ZPOT01_,
  final Array<double> D_WORK_ZPOT02_,
  final Array<double> D_WORK_ZPOT03_,
) {
  final NVAL = NVAL_.having();
  final NSVAL = NSVAL_.having();
  final A = A_.having();
  final ASAV = ASAV_.having();
  final AFAC = AFAC_.having();
  final AINV = AINV_.having();
  final B = B_.having();
  final BSAV = BSAV_.having();
  final XACT = XACT_.having();
  final X = X_.having();
  final ARF = ARF_.having();
  final ARFINV = ARFINV_.having();
  final Z_WORK_ZLATMS = Z_WORK_ZLATMS_.having();
  final Z_WORK_ZPOT02 = Z_WORK_ZPOT02_.having();
  final Z_WORK_ZPOT03 = Z_WORK_ZPOT03_.having();
  final D_WORK_ZLATMS = D_WORK_ZLATMS_.having();
  final D_WORK_ZLANHE = D_WORK_ZLANHE_.having();
  final D_WORK_ZPOT01 = D_WORK_ZPOT01_.having();
  final D_WORK_ZPOT02 = D_WORK_ZPOT02_.having();
  final D_WORK_ZPOT03 = D_WORK_ZPOT03_.having();
  const ONE = 1.0, ZERO = 0.0;
  const NTESTS = 4;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'];
  const FORMS = ['N', 'C'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  for (var IIN = 1; IIN <= NN; IIN++) {
    final N = NVAL[IIN];
    final LDA = max(N, 1);
    final LDB = max(N, 1);

    for (var IIS = 1; IIS <= NNS; IIS++) {
      final NRHS = NSVAL[IIS];

      for (var IIT = 1; IIT <= NNT; IIT++) {
        final IMAT = NTVAL[IIT];

        // If N == 0, only consider the first type

        if (N == 0 && IIT >= 1) continue;

        // Skip types 3, 4, or 5 if the matrix size is too small.

        if (IMAT == 4 && N <= 1) continue;
        if (IMAT == 5 && N <= 2) continue;

        // Do first for UPLO = 'U', then for UPLO = 'L'

        for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
          final UPLO = UPLOS[IUPLO - 1];

          // Do first for CFORM = 'N', then for CFORM = 'C'

          final RCONDC = Box(ZERO);
          for (var IFORM = 1; IFORM <= 2; IFORM++) {
            final CFORM = FORMS[IFORM - 1];

            // Set up parameters with ZLATB4 and generate a test
            // matrix with ZLATMS.

            final (TYPE: CTYPE, :KL, :KU, :ANORM, :MODE, :CNDNUM, :DIST) =
                zlatb4('ZPO', IMAT, N, N);

            srnamc.SRNAMT = 'ZLATMS';
            zlatms(N, N, DIST, ISEED, CTYPE, D_WORK_ZLATMS, MODE, CNDNUM, ANORM,
                KL, KU, UPLO, A.asMatrix(), LDA, Z_WORK_ZLATMS, INFO);

            // Check error code from ZLATMS.

            if (INFO.value != 0) {
              alaerh('ZPF', 'ZLATMS', INFO.value, 0, UPLO, N, N, -1, -1, -1,
                  IIT, NFAIL, NERRS, NOUT);
              continue;
            }

            // For types 3-5, zero one row and column of the matrix to
            // test that INFO is returned correctly.

            final ZEROT = IMAT >= 3 && IMAT <= 5;
            final int IZERO;
            if (ZEROT) {
              if (IIT == 3) {
                IZERO = 1;
              } else if (IIT == 4) {
                IZERO = N;
              } else {
                IZERO = N ~/ 2 + 1;
              }
              var IOFF = (IZERO - 1) * LDA;

              // Set row and column IZERO of A to 0.

              if (IUPLO == 1) {
                for (var I = 1; I <= IZERO - 1; I++) {
                  A[IOFF + I] = Complex.zero;
                }
                IOFF += IZERO;
                for (var I = IZERO; I <= N; I++) {
                  A[IOFF] = Complex.zero;
                  IOFF += LDA;
                }
              } else {
                IOFF = IZERO;
                for (var I = 1; I <= IZERO - 1; I++) {
                  A[IOFF] = Complex.zero;
                  IOFF += LDA;
                }
                IOFF -= IZERO;
                for (var I = IZERO; I <= N; I++) {
                  A[IOFF + I] = Complex.zero;
                }
              }
            } else {
              IZERO = 0;
            }

            // Set the imaginary part of the diagonals.

            zlaipd(N, A, LDA + 1, 0);

            // Save a copy of the matrix A in ASAV.

            zlacpy(UPLO, N, N, A.asMatrix(), LDA, ASAV.asMatrix(), LDA);

            // Compute the condition number of A (RCONDC).

            if (ZEROT) {
              RCONDC.value = ZERO;
            } else {
              // Compute the 1-norm of A.

              final ANORM =
                  zlanhe('1', UPLO, N, A.asMatrix(), LDA, D_WORK_ZLANHE);

              // Factor the matrix A.

              zpotrf(UPLO, N, A.asMatrix(), LDA, INFO);

              // Form the inverse of A.

              zpotri(UPLO, N, A.asMatrix(), LDA, INFO);

              if (N != 0) {
                // Compute the 1-norm condition number of A.

                final AINVNM =
                    zlanhe('1', UPLO, N, A.asMatrix(), LDA, D_WORK_ZLANHE);
                RCONDC.value = (ONE / ANORM) / AINVNM;

                // Restore the matrix A.

                zlacpy(UPLO, N, N, ASAV.asMatrix(), LDA, A.asMatrix(), LDA);
              }
            }

            // Form an exact solution and set the right hand side.

            srnamc.SRNAMT = 'ZLARHS';
            zlarhs('ZPO', 'N', UPLO, ' ', N, N, KL, KU, NRHS, A.asMatrix(), LDA,
                XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
            zlacpy('Full', N, NRHS, B.asMatrix(), LDA, BSAV.asMatrix(), LDA);

            // Compute the L*L' or U'*U factorization of the
            // matrix and solve the system.

            zlacpy(UPLO, N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);
            zlacpy('Full', N, NRHS, B.asMatrix(), LDB, X.asMatrix(), LDB);

            srnamc.SRNAMT = 'ZTRTTF';
            ztrttf(CFORM, UPLO, N, AFAC.asMatrix(), LDA, ARF, INFO);
            srnamc.SRNAMT = 'ZPFTRF';
            zpftrf(CFORM, UPLO, N, ARF, INFO);

            // Check error code from ZPFTRF.

            if (INFO.value != IZERO) {
              // LANGOU: there is a small hick here: IZERO should
              // always be INFO however if INFO is ZERO, ALAERH does not
              // complain.

              alaerh('ZPF', 'ZPFSV ', INFO.value, IZERO, UPLO, N, N, -1, -1,
                  NRHS, IIT, NFAIL, NERRS, NOUT);
              continue;
            }

            // Skip the tests if INFO is not 0.

            if (INFO.value != 0) {
              continue;
            }

            srnamc.SRNAMT = 'ZPFTRS';
            zpftrs(CFORM, UPLO, N, NRHS, ARF, X.asMatrix(), LDB, INFO);

            srnamc.SRNAMT = 'ZTFTTR';
            ztfttr(CFORM, UPLO, N, ARF, AFAC.asMatrix(), LDA, INFO);

            // Reconstruct matrix from factors and compute
            // residual.

            zlacpy(UPLO, N, N, AFAC.asMatrix(), LDA, ASAV.asMatrix(), LDA);
            zpot01(UPLO, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA,
                D_WORK_ZPOT01, RESULT(1));
            zlacpy(UPLO, N, N, ASAV.asMatrix(), LDA, AFAC.asMatrix(), LDA);

            // Form the inverse and compute the residual.

            if ((N % 2) == 0) {
              zlacpy('A', N + 1, N ~/ 2, ARF.asMatrix(), N + 1,
                  ARFINV.asMatrix(), N + 1);
            } else {
              zlacpy('A', N, (N + 1) ~/ 2, ARF.asMatrix(), N, ARFINV.asMatrix(),
                  N);
            }

            srnamc.SRNAMT = 'ZPFTRI';
            zpftri(CFORM, UPLO, N, ARFINV, INFO);

            srnamc.SRNAMT = 'ZTFTTR';
            ztfttr(CFORM, UPLO, N, ARFINV, AINV.asMatrix(), LDA, INFO);

            // Check error code from ZPFTRI.

            if (INFO.value != 0) {
              alaerh('ZPO', 'ZPFTRI', INFO.value, 0, UPLO, N, N, -1, -1, -1,
                  IMAT, NFAIL, NERRS, NOUT);
            }

            zpot03(
                UPLO,
                N,
                A.asMatrix(),
                LDA,
                AINV.asMatrix(),
                LDA,
                Z_WORK_ZPOT03.asMatrix(),
                LDA,
                D_WORK_ZPOT03,
                RCONDC,
                RESULT(2));

            // Compute residual of the computed solution.

            zlacpy('Full', N, NRHS, B.asMatrix(), LDA, Z_WORK_ZPOT02.asMatrix(),
                LDA);
            zpot02(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                Z_WORK_ZPOT02.asMatrix(), LDA, D_WORK_ZPOT02, RESULT(3));

            // Check solution from generated exact solution.

            zget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                RCONDC.value, RESULT(4));
            const NT = 4;

            // Print information about the tests that did not
            // pass the threshold.

            for (var K = 1; K <= NT; K++) {
              if (RESULT[K] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, 'ZPF');
                NOUT.println(
                    ' ZPFSV, UPLO=\'${UPLO.a1}\', N =${N.i5}, type ${IMAT.i1}, test(${K.i1})=${RESULT[K].g12_5}');
                NFAIL++;
              }
            }
            NRUN += NT;
          }
        }
      }
    }
  }

  // Print a summary of the results.

  alasvm('ZPF', NOUT, NFAIL, NRUN, NERRS.value);
}
