// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/zgerq2.dart';
import 'package:dart_lapack/src/zgerqf.dart';
import 'package:dart_lapack/src/zungr2.dart';
import 'package:dart_lapack/src/zungrq.dart';
import 'package:dart_lapack/src/zunmr2.dart';
import 'package:dart_lapack/src/zunmrq.dart';

import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';
import 'zgerqs.dart';

void zerrrq(final String PATH, final Nout NUNIT) {
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

  // Error exits for RQ factorization

  // ZGERQF

  srnamc.SRNAMT = 'ZGERQF';
  infoc.INFOT = 1;
  zgerqf(-1, 0, A, 1, B, W, 1, INFO);
  chkxer('ZGERQF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zgerqf(0, -1, A, 1, B, W, 1, INFO);
  chkxer('ZGERQF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  zgerqf(2, 1, A, 1, B, W, 2, INFO);
  chkxer('ZGERQF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  zgerqf(2, 1, A, 2, B, W, 1, INFO);
  chkxer('ZGERQF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // ZGERQ2

  srnamc.SRNAMT = 'ZGERQ2';
  infoc.INFOT = 1;
  zgerq2(-1, 0, A, 1, B, W, INFO);
  chkxer('ZGERQ2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zgerq2(0, -1, A, 1, B, W, INFO);
  chkxer('ZGERQ2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  zgerq2(2, 1, A, 1, B, W, INFO);
  chkxer('ZGERQ2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // ZGERQS

  srnamc.SRNAMT = 'ZGERQS';
  infoc.INFOT = 1;
  zgerqs(-1, 0, 0, A, 1, X, B.asMatrix(), 1, W, 1, INFO);
  chkxer('ZGERQS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zgerqs(0, -1, 0, A, 1, X, B.asMatrix(), 1, W, 1, INFO);
  chkxer('ZGERQS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zgerqs(2, 1, 0, A, 2, X, B.asMatrix(), 1, W, 1, INFO);
  chkxer('ZGERQS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  zgerqs(0, 0, -1, A, 1, X, B.asMatrix(), 1, W, 1, INFO);
  chkxer('ZGERQS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zgerqs(2, 2, 0, A, 1, X, B.asMatrix(), 2, W, 1, INFO);
  chkxer('ZGERQS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  zgerqs(2, 2, 0, A, 2, X, B.asMatrix(), 1, W, 1, INFO);
  chkxer('ZGERQS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 10;
  zgerqs(1, 1, 2, A, 1, X, B.asMatrix(), 1, W, 1, INFO);
  chkxer('ZGERQS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // ZUNGRQ

  srnamc.SRNAMT = 'ZUNGRQ';
  infoc.INFOT = 1;
  zungrq(-1, 0, 0, A, 1, X, W, 1, INFO);
  chkxer('ZUNGRQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zungrq(0, -1, 0, A, 1, X, W, 1, INFO);
  chkxer('ZUNGRQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zungrq(2, 1, 0, A, 2, X, W, 2, INFO);
  chkxer('ZUNGRQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  zungrq(0, 0, -1, A, 1, X, W, 1, INFO);
  chkxer('ZUNGRQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  zungrq(1, 2, 2, A, 1, X, W, 1, INFO);
  chkxer('ZUNGRQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zungrq(2, 2, 0, A, 1, X, W, 2, INFO);
  chkxer('ZUNGRQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  zungrq(2, 2, 0, A, 2, X, W, 1, INFO);
  chkxer('ZUNGRQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // ZUNGR2

  srnamc.SRNAMT = 'ZUNGR2';
  infoc.INFOT = 1;
  zungr2(-1, 0, 0, A, 1, X, W, INFO);
  chkxer('ZUNGR2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zungr2(0, -1, 0, A, 1, X, W, INFO);
  chkxer('ZUNGR2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zungr2(2, 1, 0, A, 2, X, W, INFO);
  chkxer('ZUNGR2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  zungr2(0, 0, -1, A, 1, X, W, INFO);
  chkxer('ZUNGR2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  zungr2(1, 2, 2, A, 2, X, W, INFO);
  chkxer('ZUNGR2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zungr2(2, 2, 0, A, 1, X, W, INFO);
  chkxer('ZUNGR2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // ZUNMRQ

  srnamc.SRNAMT = 'ZUNMRQ';
  infoc.INFOT = 1;
  zunmrq('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('ZUNMRQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zunmrq('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('ZUNMRQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  zunmrq('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('ZUNMRQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  zunmrq('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('ZUNMRQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zunmrq('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('ZUNMRQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zunmrq('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('ZUNMRQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zunmrq('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('ZUNMRQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  zunmrq('L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, 1, INFO);
  chkxer('ZUNMRQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  zunmrq('R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('ZUNMRQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 10;
  zunmrq('L', 'N', 2, 1, 0, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('ZUNMRQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 12;
  zunmrq('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('ZUNMRQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 12;
  zunmrq('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO);
  chkxer('ZUNMRQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // ZUNMR2

  srnamc.SRNAMT = 'ZUNMR2';
  infoc.INFOT = 1;
  zunmr2('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO);
  chkxer('ZUNMR2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zunmr2('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO);
  chkxer('ZUNMR2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  zunmr2('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO);
  chkxer('ZUNMR2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  zunmr2('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO);
  chkxer('ZUNMR2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zunmr2('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO);
  chkxer('ZUNMR2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zunmr2('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO);
  chkxer('ZUNMR2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zunmr2('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO);
  chkxer('ZUNMR2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  zunmr2('L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, INFO);
  chkxer('ZUNMR2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  zunmr2('R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, INFO);
  chkxer('ZUNMR2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 10;
  zunmr2('L', 'N', 2, 1, 0, A, 1, X, AF, 1, W, INFO);
  chkxer('ZUNMR2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // Print a summary line.

  alaesm(PATH, infoc.OK.value, NOUT);
}
