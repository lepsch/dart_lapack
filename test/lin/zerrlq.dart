// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/zgelq2.dart';
import 'package:dart_lapack/src/zgelqf.dart';
import 'package:dart_lapack/src/zungl2.dart';
import 'package:dart_lapack/src/zunglq.dart';
import 'package:dart_lapack/src/zunml2.dart';
import 'package:dart_lapack/src/zunmlq.dart';

import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void zerrlq(final String PATH, final Nout NUNIT) {
  const NMAX = 2;
  final A = Matrix<Complex>(NMAX, NMAX),
      AF = Matrix<Complex>(NMAX, NMAX),
      B = Array<Complex>(NMAX),
      W = Array<Complex>(NMAX),
      X = Array<Complex>(NMAX);
  final INFO = Box(0);

  final NOUT = NUNIT;
  NOUT.println();

  // Set the variables to innocuous values.

  for (var J = 1; J <= NMAX; J++) {
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = Complex(1.0 / (I + J), -1.0 / (I + J));
      AF[I][J] = Complex(1.0 / (I + J), -1.0 / (I + J));
    }
    B[J] = Complex.zero;
    W[J] = Complex.zero;
    X[J] = Complex.zero;
  }
  infoc.OK.value = true;

  // Error exits for LQ factorization

  // ZGELQF

  srnamc.SRNAMT = 'ZGELQF';
  infoc.INFOT = 1;
  zgelqf(-1, 0, A, 1, B, W, 1, INFO);
  chkxer('ZGELQF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zgelqf(0, -1, A, 1, B, W, 1, INFO);
  chkxer('ZGELQF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  zgelqf(2, 1, A, 1, B, W, 2, INFO);
  chkxer('ZGELQF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  zgelqf(2, 1, A, 2, B, W, 1, INFO);
  chkxer('ZGELQF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // ZGELQ2

  srnamc.SRNAMT = 'ZGELQ2';
  infoc.INFOT = 1;
  zgelq2(-1, 0, A, 1, B, W, INFO);
  chkxer('ZGELQ2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zgelq2(0, -1, A, 1, B, W, INFO);
  chkxer('ZGELQ2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  zgelq2(2, 1, A, 1, B, W, INFO);
  chkxer('ZGELQ2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // ZUNGLQ

  srnamc.SRNAMT = 'ZUNGLQ';
  infoc.INFOT = 1;
  zunglq(-1, 0, 0, A, 1, X, W, 1, INFO);
  chkxer('ZUNGLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zunglq(0, -1, 0, A, 1, X, W, 1, INFO);
  chkxer('ZUNGLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zunglq(2, 1, 0, A, 2, X, W, 2, INFO);
  chkxer('ZUNGLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  zunglq(0, 0, -1, A, 1, X, W, 1, INFO);
  chkxer('ZUNGLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  zunglq(1, 1, 2, A, 1, X, W, 1, INFO);
  chkxer('ZUNGLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zunglq(2, 2, 0, A, 1, X, W, 2, INFO);
  chkxer('ZUNGLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  zunglq(2, 2, 0, A, 2, X, W, 1, INFO);
  chkxer('ZUNGLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // ZUNGL2

  srnamc.SRNAMT = 'ZUNGL2';
  infoc.INFOT = 1;
  zungl2(-1, 0, 0, A, 1, X, W, INFO);
  chkxer('ZUNGL2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zungl2(0, -1, 0, A, 1, X, W, INFO);
  chkxer('ZUNGL2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zungl2(2, 1, 0, A, 2, X, W, INFO);
  chkxer('ZUNGL2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  zungl2(0, 0, -1, A, 1, X, W, INFO);
  chkxer('ZUNGL2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  zungl2(1, 1, 2, A, 1, X, W, INFO);
  chkxer('ZUNGL2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zungl2(2, 2, 0, A, 1, X, W, INFO);
  chkxer('ZUNGL2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // ZUNMLQ

  srnamc.SRNAMT = 'ZUNMLQ';
  infoc.INFOT = 1;
  zunmlq('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zunmlq('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  zunmlq('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  zunmlq('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zunmlq('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zunmlq('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zunmlq('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  zunmlq('L', 'N', 2, 0, 2, A, 1, X, AF, 2, W, 1, INFO);
  chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  zunmlq('R', 'N', 0, 2, 2, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 10;
  zunmlq('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO);
  chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 12;
  zunmlq('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 12;
  zunmlq('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO);
  chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // ZUNML2

  srnamc.SRNAMT = 'ZUNML2';
  infoc.INFOT = 1;
  zunml2('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO);
  chkxer('ZUNML2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zunml2('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO);
  chkxer('ZUNML2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  zunml2('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO);
  chkxer('ZUNML2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  zunml2('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO);
  chkxer('ZUNML2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zunml2('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO);
  chkxer('ZUNML2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zunml2('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO);
  chkxer('ZUNML2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zunml2('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO);
  chkxer('ZUNML2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  zunml2('L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, INFO);
  chkxer('ZUNML2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  zunml2('R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, INFO);
  chkxer('ZUNML2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 10;
  zunml2('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO);
  chkxer('ZUNML2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // Print a summary line.

  alaesm(PATH, infoc.OK.value, NOUT);
}
