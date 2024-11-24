// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/lapack.dart';

import '../test_driver.dart';
import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void derrlqt(final String PATH, final Nout NUNIT, final TestDriver test) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 2;
  final A = Matrix<double>(NMAX, NMAX),
      T = Matrix<double>(NMAX, NMAX),
      W = Array<double>(NMAX),
      C = Matrix<double>(NMAX, NMAX);
  final INFO = Box(0);

  final NOUT = infoc.NOUT = NUNIT;
  final OK = infoc.OK;
  final LERR = infoc.LERR;
  infoc.NOUT.println();

  // Set the variables to innocuous values.

  for (var J = 1; J <= NMAX; J++) {
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = 1.0 / (I + J);
      C[I][J] = 1.0 / (I + J);
      T[I][J] = 1.0 / (I + J);
    }
    W[J] = 0.0;
  }
  infoc.OK.value = true;

  // Error exits for LQT factorization

  test('DGELQT', () {
    srnamc.SRNAMT = 'DGELQT';
    infoc.INFOT = 1;
    dgelqt(-1, 0, 1, A, 1, T, 1, W, INFO);
    chkxer('DGELQT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 2;
    dgelqt(0, -1, 1, A, 1, T, 1, W, INFO);
    chkxer('DGELQT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 3;
    dgelqt(0, 0, 0, A, 1, T, 1, W, INFO);
    chkxer('DGELQT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 5;
    dgelqt(2, 1, 1, A, 1, T, 1, W, INFO);
    chkxer('DGELQT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 7;
    dgelqt(2, 2, 2, A, 2, T, 1, W, INFO);
    chkxer('DGELQT', infoc.INFOT, NOUT, LERR, OK, test);
  });

  test('DGELQT3', () {
    srnamc.SRNAMT = 'DGELQT3';
    infoc.INFOT = 1;
    dgelqt3(-1, 0, A, 1, T, 1, INFO);
    chkxer('DGELQT3', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 2;
    dgelqt3(0, -1, A, 1, T, 1, INFO);
    chkxer('DGELQT3', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 4;
    dgelqt3(2, 2, A, 1, T, 1, INFO);
    chkxer('DGELQT3', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 6;
    dgelqt3(2, 2, A, 2, T, 1, INFO);
    chkxer('DGELQT3', infoc.INFOT, NOUT, LERR, OK, test);
  });

  test('DGEMLQT', () {
    srnamc.SRNAMT = 'DGEMLQT';
    infoc.INFOT = 1;
    dgemlqt('/', 'N', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO);
    chkxer('DGEMLQT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 2;
    dgemlqt('L', '/', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO);
    chkxer('DGEMLQT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 3;
    dgemlqt('L', 'N', -1, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO);
    chkxer('DGEMLQT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 4;
    dgemlqt('L', 'N', 0, -1, 0, 1, A, 1, T, 1, C, 1, W, INFO);
    chkxer('DGEMLQT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 5;
    dgemlqt('L', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO);
    chkxer('DGEMLQT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 5;
    dgemlqt('R', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO);
    chkxer('DGEMLQT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 6;
    dgemlqt('L', 'N', 0, 0, 0, 0, A, 1, T, 1, C, 1, W, INFO);
    chkxer('DGEMLQT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 8;
    dgemlqt('R', 'N', 2, 2, 2, 1, A, 1, T, 1, C, 1, W, INFO);
    chkxer('DGEMLQT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 8;
    dgemlqt('L', 'N', 2, 2, 2, 1, A, 1, T, 1, C, 1, W, INFO);
    chkxer('DGEMLQT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 10;
    dgemlqt('R', 'N', 1, 1, 1, 1, A, 1, T, 0, C, 1, W, INFO);
    chkxer('DGEMLQT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 12;
    dgemlqt('L', 'N', 1, 1, 1, 1, A, 1, T, 1, C, 0, W, INFO);
    chkxer('DGEMLQT', infoc.INFOT, NOUT, LERR, OK, test);
  });

  // Print a summary line.
  alaesm(PATH, OK.value, NOUT);
}
