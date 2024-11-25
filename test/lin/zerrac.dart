// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/zcposv.dart';

import 'chkxer.dart';
import 'common.dart';

void zerrac(Nout NUNIT) {
  const NMAX = 4;
  final A = Matrix<Complex>(NMAX, NMAX),
      AF = Matrix<Complex>(NMAX, NMAX),
      B = Array<Complex>(NMAX),
      C = Array<Complex>(NMAX),
      R = Array<Complex>(NMAX),
      R1 = Array<Complex>(NMAX),
      R2 = Array<Complex>(NMAX),
      W = Array<Complex>(2 * NMAX),
      X = Array<Complex>(NMAX);
  final RWORK = Array<double>(NMAX);
  final WORK = Array<Complex>(NMAX * NMAX);
  final SWORK = Array<Complex>(NMAX * NMAX);

  final NOUT = NUNIT;
  NOUT.println();

  // Set the variables to innocuous values.

  for (var J = 1; J <= NMAX; J++) {
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = (1.0 / (I + J)).toComplex();
      AF[I][J] = (1.0 / (I + J)).toComplex();
    }
    B[J] = Complex.zero;
    R1[J] = Complex.zero;
    R2[J] = Complex.zero;
    W[J] = Complex.zero;
    X[J] = Complex.zero;
    C[J] = Complex.zero;
    R[J] = Complex.zero;
  }
  infoc.OK.value = true;

  final INFO = Box(0), ITER = Box(0);
  srnamc.SRNAMT = 'ZCPOSV';
  infoc.INFOT = 1;
  zcposv('/', 0, 0, A, 1, B.asMatrix(), 1, X.asMatrix(), 1, WORK.asMatrix(),
      SWORK, RWORK, ITER, INFO);
  chkxer('ZCPOSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zcposv('U', -1, 0, A, 1, B.asMatrix(), 1, X.asMatrix(), 1, WORK.asMatrix(),
      SWORK, RWORK, ITER, INFO);
  chkxer('ZCPOSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  zcposv('U', 0, -1, A, 1, B.asMatrix(), 1, X.asMatrix(), 1, WORK.asMatrix(),
      SWORK, RWORK, ITER, INFO);
  chkxer('ZCPOSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zcposv('U', 2, 1, A, 1, B.asMatrix(), 2, X.asMatrix(), 2, WORK.asMatrix(),
      SWORK, RWORK, ITER, INFO);
  chkxer('ZCPOSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  zcposv('U', 2, 1, A, 2, B.asMatrix(), 1, X.asMatrix(), 2, WORK.asMatrix(),
      SWORK, RWORK, ITER, INFO);
  chkxer('ZCPOSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 9;
  zcposv('U', 2, 1, A, 2, B.asMatrix(), 2, X.asMatrix(), 1, WORK.asMatrix(),
      SWORK, RWORK, ITER, INFO);
  chkxer('ZCPOSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // Print a summary line.

  if (infoc.OK.value) {
    NOUT.println(' ZCPOSV drivers passed the tests of the error exits');
  } else {
    NOUT.println(' *** ZCPOSV drivers failed the tests of the error exits ***');
  }
}
