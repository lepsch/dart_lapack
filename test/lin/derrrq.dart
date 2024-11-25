// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';

import '../test_driver.dart';
import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';
import 'dgerqs.dart';

void derrrq(final String PATH, final Nout NUNIT, final TestDriver test) {
  const NMAX = 2;
  final A = Matrix<double>(NMAX, NMAX),
      AF = Matrix<double>(NMAX, NMAX),
      B = Array<double>(NMAX),
      W = Array<double>(NMAX),
      X = Array<double>(NMAX);
  final INFO = Box(0);

  final NOUT = infoc.NOUT = NUNIT;
  final OK = infoc.OK;
  final LERR = infoc.LERR;
  NOUT.println();

  // Set the variables to innocuous values.

  for (var J = 1; J <= NMAX; J++) {
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = 1.0 / (I + J);
      AF[I][J] = 1.0 / (I + J);
    }
    B[J] = 0.0;
    W[J] = 0.0;
    X[J] = 0.0;
  }
  OK.value = true;

  // Error exits for RQ factorization

  test('DGERQF', () {
    srnamc.SRNAMT = 'DGERQF';
    infoc.INFOT = 1;
    dgerqf(-1, 0, A, 1, B, W, 1, INFO);
    chkxer('DGERQF', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dgerqf(0, -1, A, 1, B, W, 1, INFO);
    chkxer('DGERQF', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 4;
    dgerqf(2, 1, A, 1, B, W, 2, INFO);
    chkxer('DGERQF', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 7;
    dgerqf(2, 1, A, 2, B, W, 1, INFO);
    chkxer('DGERQF', infoc.INFOT, NOUT, LERR, OK);
  });

  test('DGERQ2', () {
    srnamc.SRNAMT = 'DGERQ2';
    infoc.INFOT = 1;
    dgerq2(-1, 0, A, 1, B, W, INFO);
    chkxer('DGERQ2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dgerq2(0, -1, A, 1, B, W, INFO);
    chkxer('DGERQ2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 4;
    dgerq2(2, 1, A, 1, B, W, INFO);
    chkxer('DGERQ2', infoc.INFOT, NOUT, LERR, OK);
  });

  test('DGERQS', () {
    srnamc.SRNAMT = 'DGERQS';
    infoc.INFOT = 1;
    dgerqs(-1, 0, 0, A, 1, X, B.asMatrix(), 1, W, 1, INFO);
    chkxer('DGERQS', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dgerqs(0, -1, 0, A, 1, X, B.asMatrix(), 1, W, 1, INFO);
    chkxer('DGERQS', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dgerqs(2, 1, 0, A, 2, X, B.asMatrix(), 1, W, 1, INFO);
    chkxer('DGERQS', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 3;
    dgerqs(0, 0, -1, A, 1, X, B.asMatrix(), 1, W, 1, INFO);
    chkxer('DGERQS', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dgerqs(2, 2, 0, A, 1, X, B.asMatrix(), 2, W, 1, INFO);
    chkxer('DGERQS', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 8;
    dgerqs(2, 2, 0, A, 2, X, B.asMatrix(), 1, W, 1, INFO);
    chkxer('DGERQS', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 10;
    dgerqs(1, 1, 2, A, 1, X, B.asMatrix(), 1, W, 1, INFO);
    chkxer('DGERQS', infoc.INFOT, NOUT, LERR, OK);
  });

  test('DORGRQ', () {
    srnamc.SRNAMT = 'DORGRQ';
    infoc.INFOT = 1;
    dorgrq(-1, 0, 0, A, 1, X, W, 1, INFO);
    chkxer('DORGRQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dorgrq(0, -1, 0, A, 1, X, W, 1, INFO);
    chkxer('DORGRQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dorgrq(2, 1, 0, A, 2, X, W, 2, INFO);
    chkxer('DORGRQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 3;
    dorgrq(0, 0, -1, A, 1, X, W, 1, INFO);
    chkxer('DORGRQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 3;
    dorgrq(1, 2, 2, A, 1, X, W, 1, INFO);
    chkxer('DORGRQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dorgrq(2, 2, 0, A, 1, X, W, 2, INFO);
    chkxer('DORGRQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 8;
    dorgrq(2, 2, 0, A, 2, X, W, 1, INFO);
    chkxer('DORGRQ', infoc.INFOT, NOUT, LERR, OK);
  });

  test('DORGR2', () {
    srnamc.SRNAMT = 'DORGR2';
    infoc.INFOT = 1;
    dorgr2(-1, 0, 0, A, 1, X, W, INFO);
    chkxer('DORGR2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dorgr2(0, -1, 0, A, 1, X, W, INFO);
    chkxer('DORGR2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dorgr2(2, 1, 0, A, 2, X, W, INFO);
    chkxer('DORGR2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 3;
    dorgr2(0, 0, -1, A, 1, X, W, INFO);
    chkxer('DORGR2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 3;
    dorgr2(1, 2, 2, A, 2, X, W, INFO);
    chkxer('DORGR2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dorgr2(2, 2, 0, A, 1, X, W, INFO);
    chkxer('DORGR2', infoc.INFOT, NOUT, LERR, OK);
  });

  test('DORMRQ', () {
    srnamc.SRNAMT = 'DORMRQ';
    infoc.INFOT = 1;
    dormrq('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMRQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dormrq('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMRQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 3;
    dormrq('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMRQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 4;
    dormrq('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMRQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dormrq('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMRQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dormrq('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMRQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dormrq('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMRQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 7;
    dormrq('L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, 1, INFO);
    chkxer('DORMRQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 7;
    dormrq('R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMRQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 10;
    dormrq('L', 'N', 2, 1, 0, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMRQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 12;
    dormrq('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMRQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 12;
    dormrq('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO);
    chkxer('DORMRQ', infoc.INFOT, NOUT, LERR, OK);
  });

  test('DORMR2', () {
    srnamc.SRNAMT = 'DORMR2';
    infoc.INFOT = 1;
    dormr2('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO);
    chkxer('DORMR2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dormr2('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO);
    chkxer('DORMR2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 3;
    dormr2('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO);
    chkxer('DORMR2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 4;
    dormr2('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO);
    chkxer('DORMR2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dormr2('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO);
    chkxer('DORMR2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dormr2('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO);
    chkxer('DORMR2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dormr2('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO);
    chkxer('DORMR2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 7;
    dormr2('L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, INFO);
    chkxer('DORMR2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 7;
    dormr2('R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, INFO);
    chkxer('DORMR2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 10;
    dormr2('L', 'N', 2, 1, 0, A, 1, X, AF, 1, W, INFO);
    chkxer('DORMR2', infoc.INFOT, NOUT, LERR, OK);
  });

  // Print a summary line.
  alaesm(PATH, infoc.OK.value, infoc.NOUT);
}
