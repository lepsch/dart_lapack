// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/zlanhe.dart';
import 'package:dart_lapack/src/zlanhf.dart';
import 'package:dart_lapack/src/ztrttf.dart';

import '../matgen/zlarnd.dart';
import 'common.dart';

void zdrvrf1(
  final Nout NOUT,
  final int NN,
  final Array<int> NVAL_,
  final double THRESH,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> ARF_,
  final Array<double> WORK_,
) {
  final NVAL = NVAL_.having();
  final ARF = ARF_.having();
  final A = A_.having(ld: LDA);
  final WORK = WORK_.having();
  const ONE = 1.0;
  const NTESTS = 1;
  final RESULT = Array<double>(NTESTS);
  final ISEED = Array<int>(4);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'];
  const FORMS = ['N', 'C'];
  const NORMS = ['M', '1', 'I', 'F'];

  // Initialize constants and the random number seed.

  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  final INFO = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  final EPS = dlamch('Precision');
  var SMALL = dlamch('Safe minimum');
  var LARGE = ONE / SMALL;
  SMALL *= LDA * LDA;
  LARGE /= LDA / LDA;

  for (var IIN = 1; IIN <= NN; IIN++) {
    final N = NVAL[IIN];

    for (var IIT = 1; IIT <= 3; IIT++) {
      // Nothing to do for N=0
      if (N == 0) break;

      // IIT = 1 : random matrix
      // IIT = 2 : random matrix scaled near underflow
      // IIT = 3 : random matrix scaled near overflow

      for (var J = 1; J <= N; J++) {
        for (var I = 1; I <= N; I++) {
          A[I][J] = zlarnd(4, ISEED);
        }
      }

      if (IIT == 2) {
        for (var J = 1; J <= N; J++) {
          for (var I = 1; I <= N; I++) {
            A[I][J] *= LARGE.toComplex();
          }
        }
      }

      if (IIT == 3) {
        for (var J = 1; J <= N; J++) {
          for (var I = 1; I <= N; I++) {
            A[I][J] *= SMALL.toComplex();
          }
        }
      }

      // Do first for UPLO = 'U', then for UPLO = 'L'

      for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
        final UPLO = UPLOS[IUPLO - 1];

        // Do first for CFORM = 'N', then for CFORM = 'C'

        for (var IFORM = 1; IFORM <= 2; IFORM++) {
          final CFORM = FORMS[IFORM - 1];

          srnamc.SRNAMT = 'ZTRTTF';
          ztrttf(CFORM, UPLO, N, A, LDA, ARF, INFO);

          // Check error code from ZTRTTF

          if (INFO.value != 0) {
            if (NFAIL == 0 && NERRS.value == 0) {
              NOUT.println();
              NOUT.print9999();
            }
            NOUT.println(
                '      Error in ${srnamc.SRNAMT.a6} with UPLO=\'${UPLO.a1}\', FORM=\'${CFORM.a1}\', N=${N.i5}');
            NERRS.value++;
            continue;
          }

          for (var INORM = 1; INORM <= 4; INORM++) {
            // Check all four norms: 'M', '1', 'I', 'F'

            final NORM = NORMS[INORM];
            final NORMARF = zlanhf(NORM, CFORM, UPLO, N, ARF, WORK);
            final NORMA = zlanhe(NORM, UPLO, N, A, LDA, WORK);

            RESULT[1] = (NORMA - NORMARF) / NORMA / EPS;
            NRUN++;

            if (RESULT[1] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) {
                NOUT.println();
                NOUT.print9999();
              }
              NOUT.println(
                  '      Failure in ZLANHF N=${N.i5} TYPE=${IIT.i5} UPLO=\'${UPLO.a1}\', FORM =\'${CFORM.a1}\', NORM=\'${NORM.a1}\', test=${RESULT[1].g12_5}');
              NFAIL++;
            }
          }
        }
      }
    }
  }

  // Print a summary of the results.

  if (NFAIL == 0) {
    NOUT.println(
        ' All tests for ZLANHF auxiliary routine passed the threshold ( ${NRUN.i5} tests run)');
  } else {
    NOUT.println(
        ' ZLANHF auxiliary routine:${NFAIL.i5} out of ${NRUN.i5} tests failed to pass the threshold');
  }
  if (NERRS.value != 0) {
    NOUT.println(
        '${' ' * 26}${NERRS.value.i5} error message recorded (ZLANHF)');
  }
}

extension on Nout {
  void print9999() {
    println('  *** Error(s) or Failure(s) while testing ZLANHF ***');
  }
}
