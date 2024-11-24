// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/lapack.dart';

import '../eig/chkxer.dart';
import '../test_driver.dart';
import 'alaesm.dart';
import 'common.dart';
import 'dgeqls.dart';

void derrql(final String PATH, final Nout NUNIT, final TestDriver test) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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

  // Error exits for QL factorization

  test('DGEQLF', () {
    srnamc.SRNAMT = 'DGEQLF';
    infoc.INFOT = 1;
    dgeqlf(-1, 0, A, 1, B, W, 1, INFO);
    chkxer('DGEQLF', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dgeqlf(0, -1, A, 1, B, W, 1, INFO);
    chkxer('DGEQLF', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 4;
    dgeqlf(2, 1, A, 1, B, W, 1, INFO);
    chkxer('DGEQLF', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 7;
    dgeqlf(1, 2, A, 1, B, W, 1, INFO);
    chkxer('DGEQLF', infoc.INFOT, NOUT, LERR, OK);
  });

  test('DGEQL2', () {
    srnamc.SRNAMT = 'DGEQL2';
    infoc.INFOT = 1;
    dgeql2(-1, 0, A, 1, B, W, INFO);
    chkxer('DGEQL2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dgeql2(0, -1, A, 1, B, W, INFO);
    chkxer('DGEQL2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 4;
    dgeql2(2, 1, A, 1, B, W, INFO);
    chkxer('DGEQL2', infoc.INFOT, NOUT, LERR, OK);
  });

  test('DGEQLS', () {
    srnamc.SRNAMT = 'DGEQLS';
    infoc.INFOT = 1;
    dgeqls(-1, 0, 0, A, 1, X, B.asMatrix(), 1, W, 1, INFO);
    chkxer('DGEQLS', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dgeqls(0, -1, 0, A, 1, X, B.asMatrix(), 1, W, 1, INFO);
    chkxer('DGEQLS', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dgeqls(1, 2, 0, A, 1, X, B.asMatrix(), 1, W, 1, INFO);
    chkxer('DGEQLS', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 3;
    dgeqls(0, 0, -1, A, 1, X, B.asMatrix(), 1, W, 1, INFO);
    chkxer('DGEQLS', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dgeqls(2, 1, 0, A, 1, X, B.asMatrix(), 2, W, 1, INFO);
    chkxer('DGEQLS', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 8;
    dgeqls(2, 1, 0, A, 2, X, B.asMatrix(), 1, W, 1, INFO);
    chkxer('DGEQLS', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 10;
    dgeqls(1, 1, 2, A, 1, X, B.asMatrix(), 1, W, 1, INFO);
    chkxer('DGEQLS', infoc.INFOT, NOUT, LERR, OK);
  });

  test('DORGQL', () {
    srnamc.SRNAMT = 'DORGQL';
    infoc.INFOT = 1;
    dorgql(-1, 0, 0, A, 1, X, W, 1, INFO);
    chkxer('DORGQL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dorgql(0, -1, 0, A, 1, X, W, 1, INFO);
    chkxer('DORGQL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dorgql(1, 2, 0, A, 1, X, W, 2, INFO);
    chkxer('DORGQL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 3;
    dorgql(0, 0, -1, A, 1, X, W, 1, INFO);
    chkxer('DORGQL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 3;
    dorgql(1, 1, 2, A, 1, X, W, 1, INFO);
    chkxer('DORGQL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dorgql(2, 1, 0, A, 1, X, W, 1, INFO);
    chkxer('DORGQL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 8;
    dorgql(2, 2, 0, A, 2, X, W, 1, INFO);
    chkxer('DORGQL', infoc.INFOT, NOUT, LERR, OK);
  });

  test('DORG2L', () {
    srnamc.SRNAMT = 'DORG2L';
    infoc.INFOT = 1;
    dorg2l(-1, 0, 0, A, 1, X, W, INFO);
    chkxer('DORG2L', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dorg2l(0, -1, 0, A, 1, X, W, INFO);
    chkxer('DORG2L', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dorg2l(1, 2, 0, A, 1, X, W, INFO);
    chkxer('DORG2L', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 3;
    dorg2l(0, 0, -1, A, 1, X, W, INFO);
    chkxer('DORG2L', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 3;
    dorg2l(2, 1, 2, A, 2, X, W, INFO);
    chkxer('DORG2L', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dorg2l(2, 1, 0, A, 1, X, W, INFO);
    chkxer('DORG2L', infoc.INFOT, NOUT, LERR, OK);
  });

  test('DORMQL', () {
    srnamc.SRNAMT = 'DORMQL';
    infoc.INFOT = 1;
    dormql('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMQL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dormql('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMQL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 3;
    dormql('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMQL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 4;
    dormql('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMQL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dormql('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMQL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dormql('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMQL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dormql('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMQL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 7;
    dormql('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO);
    chkxer('DORMQL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 7;
    dormql('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMQL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 10;
    dormql('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO);
    chkxer('DORMQL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 12;
    dormql('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMQL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 12;
    dormql('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO);
    chkxer('DORMQL', infoc.INFOT, NOUT, LERR, OK);
  });

  test('DORM2L', () {
    srnamc.SRNAMT = 'DORM2L';
    infoc.INFOT = 1;
    dorm2l('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO);
    chkxer('DORM2L', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dorm2l('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO);
    chkxer('DORM2L', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 3;
    dorm2l('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO);
    chkxer('DORM2L', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 4;
    dorm2l('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO);
    chkxer('DORM2L', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dorm2l('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO);
    chkxer('DORM2L', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dorm2l('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO);
    chkxer('DORM2L', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dorm2l('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO);
    chkxer('DORM2L', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 7;
    dorm2l('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, INFO);
    chkxer('DORM2L', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 7;
    dorm2l('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, INFO);
    chkxer('DORM2L', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 10;
    dorm2l('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO);
    chkxer('DORM2L', infoc.INFOT, NOUT, LERR, OK);
  });

  // Print a summary line.
  alaesm(PATH, OK.value, NOUT);
}
