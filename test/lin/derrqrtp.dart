// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';

import '../test_driver.dart';
import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void derrqrtp(final String PATH, final Nout NUNIT, final TestDriver test) {
  const NMAX = 2;
  final A = Matrix<double>(NMAX, NMAX),
      T = Matrix<double>(NMAX, NMAX),
      W = Array<double>(NMAX),
      B = Matrix<double>(NMAX, NMAX),
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

  // Error exits for TPQRT factorization

  test('DTPQRT', () {
    srnamc.SRNAMT = 'DTPQRT';
    infoc.INFOT = 1;
    dtpqrt(-1, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO);
    chkxer('DTPQRT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 2;
    dtpqrt(1, -1, 0, 1, A, 1, B, 1, T, 1, W, INFO);
    chkxer('DTPQRT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 3;
    dtpqrt(0, 1, -1, 1, A, 1, B, 1, T, 1, W, INFO);
    chkxer('DTPQRT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 3;
    dtpqrt(0, 1, 1, 1, A, 1, B, 1, T, 1, W, INFO);
    chkxer('DTPQRT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 4;
    dtpqrt(0, 1, 0, 0, A, 1, B, 1, T, 1, W, INFO);
    chkxer('DTPQRT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 4;
    dtpqrt(0, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO);
    chkxer('DTPQRT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 6;
    dtpqrt(1, 2, 0, 2, A, 1, B, 1, T, 1, W, INFO);
    chkxer('DTPQRT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 8;
    dtpqrt(2, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO);
    chkxer('DTPQRT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 10;
    dtpqrt(2, 2, 1, 2, A, 2, B, 2, T, 1, W, INFO);
    chkxer('DTPQRT', infoc.INFOT, NOUT, LERR, OK, test);
  });

  test('DTPQRT2', () {
    srnamc.SRNAMT = 'DTPQRT2';
    infoc.INFOT = 1;
    dtpqrt2(-1, 0, 0, A, 1, B, 1, T, 1, INFO);
    chkxer('DTPQRT2', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 2;
    dtpqrt2(0, -1, 0, A, 1, B, 1, T, 1, INFO);
    chkxer('DTPQRT2', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 3;
    dtpqrt2(0, 0, -1, A, 1, B, 1, T, 1, INFO);
    chkxer('DTPQRT2', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 5;
    dtpqrt2(2, 2, 0, A, 1, B, 2, T, 2, INFO);
    chkxer('DTPQRT2', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 7;
    dtpqrt2(2, 2, 0, A, 2, B, 1, T, 2, INFO);
    chkxer('DTPQRT2', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 9;
    dtpqrt2(2, 2, 0, A, 2, B, 2, T, 1, INFO);
    chkxer('DTPQRT2', infoc.INFOT, NOUT, LERR, OK, test);
  });

  test('DTPMQRT', () {
    srnamc.SRNAMT = 'DTPMQRT';
    infoc.INFOT = 1;
    dtpmqrt('/', 'N', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
    chkxer('DTPMQRT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 2;
    dtpmqrt('L', '/', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
    chkxer('DTPMQRT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 3;
    dtpmqrt('L', 'N', -1, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
    chkxer('DTPMQRT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 4;
    dtpmqrt('L', 'N', 0, -1, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
    chkxer('DTPMQRT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 5;
    dtpmqrt('L', 'N', 0, 0, -1, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
    infoc.INFOT = 6;
    dtpmqrt('L', 'N', 0, 0, 0, -1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
    chkxer('DTPMQRT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 7;
    dtpmqrt('L', 'N', 0, 0, 0, 0, 0, A, 1, T, 1, B, 1, C, 1, W, INFO);
    chkxer('DTPMQRT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 9;
    dtpmqrt('R', 'N', 1, 2, 1, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
    chkxer('DTPMQRT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 9;
    dtpmqrt('L', 'N', 2, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
    chkxer('DTPMQRT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 11;
    dtpmqrt('R', 'N', 1, 1, 1, 1, 1, A, 1, T, 0, B, 1, C, 1, W, INFO);
    chkxer('DTPMQRT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 13;
    dtpmqrt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 0, C, 1, W, INFO);
    chkxer('DTPMQRT', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 15;
    dtpmqrt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 0, W, INFO);
    chkxer('DTPMQRT', infoc.INFOT, NOUT, LERR, OK, test);
  });

  // Print a summary line.
  alaesm(PATH, OK.value, NOUT);
}
