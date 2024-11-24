// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';

import '../test_driver.dart';
import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void derrtsqr(final String PATH, final Nout NUNIT, final TestDriver test) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 2;
  final A = Matrix<double>(NMAX, NMAX),
      T = Matrix<double>(NMAX, NMAX),
      W = Array<double>(NMAX),
      C = Matrix<double>(NMAX, NMAX),
      TAU = Array<double>(NMAX * 2);
  final INFO = Box(0);

  final NOUT = infoc.NOUT = NUNIT;
  final OK = infoc.OK;
  final LERR = infoc.LERR;
  NOUT.println();

  // Set the variables to innocuous values.

  for (var J = 1; J <= NMAX; J++) {
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = 1.0 / (I + J);
      C[I][J] = 1.0 / (I + J);
      T[I][J] = 1.0 / (I + J);
    }
    W[J] = 0.0;
  }
  OK.value = true;

  // Error exits for TS factorization

  test('DGEQR', () {
    srnamc.SRNAMT = 'DGEQR';
    infoc.INFOT = 1;
    dgeqr(-1, 0, A, 1, TAU, 1, W, 1, INFO);
    chkxer('DGEQR', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 2;
    dgeqr(0, -1, A, 1, TAU, 1, W, 1, INFO);
    chkxer('DGEQR', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 4;
    dgeqr(1, 1, A, 0, TAU, 1, W, 1, INFO);
    chkxer('DGEQR', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 6;
    dgeqr(3, 2, A, 3, TAU, 1, W, 1, INFO);
    chkxer('DGEQR', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 8;
    dgeqr(3, 2, A, 3, TAU, 7, W, 0, INFO);
    chkxer('DGEQR', infoc.INFOT, NOUT, LERR, OK, test);
  });

  test('DLATSQR', () {
    final MB = 1;
    final NB = 1;
    srnamc.SRNAMT = 'DLATSQR';
    infoc.INFOT = 1;
    dlatsqr(-1, 0, MB, NB, A, 1, TAU.asMatrix(), 1, W, 1, INFO);
    chkxer('DLATSQR', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 2;
    dlatsqr(1, 2, MB, NB, A, 1, TAU.asMatrix(), 1, W, 1, INFO);
    chkxer('DLATSQR', infoc.INFOT, NOUT, LERR, OK, test);
    dlatsqr(0, -1, MB, NB, A, 1, TAU.asMatrix(), 1, W, 1, INFO);
    chkxer('DLATSQR', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 3;
    dlatsqr(2, 1, -1, NB, A, 2, TAU.asMatrix(), 1, W, 1, INFO);
    chkxer('DLATSQR', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 4;
    dlatsqr(2, 1, MB, 2, A, 2, TAU.asMatrix(), 1, W, 1, INFO);
    chkxer('DLATSQR', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 6;
    dlatsqr(2, 1, MB, NB, A, 1, TAU.asMatrix(), 1, W, 1, INFO);
    chkxer('DLATSQR', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 8;
    dlatsqr(2, 1, MB, NB, A, 2, TAU.asMatrix(), 0, W, 1, INFO);
    chkxer('DLATSQR', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 10;
    dlatsqr(2, 1, MB, NB, A, 2, TAU.asMatrix(), 2, W, 0, INFO);
    chkxer('DLATSQR', infoc.INFOT, NOUT, LERR, OK, test);
  });

  test('DGEMQR', () {
    TAU[1] = 1;
    TAU[2] = 1;
    TAU[3] = 1;
    TAU[4] = 1;
    srnamc.SRNAMT = 'DGEMQR';
    infoc.INFOT = 1;
    dgemqr('/', 'N', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1, INFO);
    chkxer('DGEMQR', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 2;
    dgemqr('L', '/', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1, INFO);
    chkxer('DGEMQR', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 3;
    dgemqr('L', 'N', -1, 0, 0, A, 1, TAU, 1, C, 1, W, 1, INFO);
    chkxer('DGEMQR', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 4;
    dgemqr('L', 'N', 0, -1, 0, A, 1, TAU, 1, C, 1, W, 1, INFO);
    chkxer('DGEMQR', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 5;
    dgemqr('L', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W, 1, INFO);
    chkxer('DGEMQR', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 5;
    dgemqr('R', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W, 1, INFO);
    chkxer('DGEMQR', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 7;
    dgemqr('L', 'N', 2, 1, 0, A, 0, TAU, 1, C, 1, W, 1, INFO);
    chkxer('DGEMQR', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 9;
    dgemqr('R', 'N', 2, 2, 1, A, 2, TAU, 0, C, 1, W, 1, INFO);
    chkxer('DGEMQR', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 9;
    dgemqr('L', 'N', 2, 2, 1, A, 2, TAU, 0, C, 1, W, 1, INFO);
    chkxer('DGEMQR', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 11;
    dgemqr('L', 'N', 2, 1, 1, A, 2, TAU, 6, C, 0, W, 1, INFO);
    chkxer('DGEMQR', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 13;
    dgemqr('L', 'N', 2, 2, 1, A, 2, TAU, 6, C, 2, W, 0, INFO);
    chkxer('DGEMQR', infoc.INFOT, NOUT, LERR, OK, test);
  });

  test('DGELQ', () {
    srnamc.SRNAMT = 'DGELQ';
    infoc.INFOT = 1;
    dgelq(-1, 0, A, 1, TAU, 1, W, 1, INFO);
    chkxer('DGELQ', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 2;
    dgelq(0, -1, A, 1, TAU, 1, W, 1, INFO);
    chkxer('DGELQ', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 4;
    dgelq(1, 1, A, 0, TAU, 1, W, 1, INFO);
    chkxer('DGELQ', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 6;
    dgelq(2, 3, A, 3, TAU, 1, W, 1, INFO);
    chkxer('DGELQ', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 8;
    dgelq(2, 3, A, 3, TAU, 7, W, 0, INFO);
    chkxer('DGELQ', infoc.INFOT, NOUT, LERR, OK, test);
  });

  test('DLASWLQ', () {
    final MB = 1;
    final NB = 1;
    srnamc.SRNAMT = 'DLASWLQ';
    infoc.INFOT = 1;
    dlaswlq(-1, 0, MB, NB, A, 1, TAU.asMatrix(), 1, W, 1, INFO);
    chkxer('DLASWLQ', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 2;
    dlaswlq(2, 1, MB, NB, A, 1, TAU.asMatrix(), 1, W, 1, INFO);
    chkxer('DLASWLQ', infoc.INFOT, NOUT, LERR, OK, test);
    dlaswlq(0, -1, MB, NB, A, 1, TAU.asMatrix(), 1, W, 1, INFO);
    chkxer('DLASWLQ', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 3;
    dlaswlq(1, 2, -1, NB, A, 1, TAU.asMatrix(), 1, W, 1, INFO);
    chkxer('DLASWLQ', infoc.INFOT, NOUT, LERR, OK, test);
    dlaswlq(1, 1, 2, NB, A, 1, TAU.asMatrix(), 1, W, 1, INFO);
    chkxer('DLASWLQ', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 4;
    dlaswlq(1, 2, MB, -1, A, 1, TAU.asMatrix(), 1, W, 1, INFO);
    chkxer('DLASWLQ', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 6;
    dlaswlq(1, 2, MB, NB, A, 0, TAU.asMatrix(), 1, W, 1, INFO);
    chkxer('DLASWLQ', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 8;
    dlaswlq(1, 2, MB, NB, A, 1, TAU.asMatrix(), 0, W, 1, INFO);
    chkxer('DLASWLQ', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 10;
    dlaswlq(1, 2, MB, NB, A, 1, TAU.asMatrix(), 1, W, 0, INFO);
    chkxer('DLASWLQ', infoc.INFOT, NOUT, LERR, OK, test);
  });

  test('DGEMLQ', () {
    TAU[1] = 1;
    TAU[2] = 1;
    srnamc.SRNAMT = 'DGEMLQ';
    infoc.INFOT = 1;
    dgemlq('/', 'N', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1, INFO);
    chkxer('DGEMLQ', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 2;
    dgemlq('L', '/', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1, INFO);
    chkxer('DGEMLQ', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 3;
    dgemlq('L', 'N', -1, 0, 0, A, 1, TAU, 1, C, 1, W, 1, INFO);
    chkxer('DGEMLQ', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 4;
    dgemlq('L', 'N', 0, -1, 0, A, 1, TAU, 1, C, 1, W, 1, INFO);
    chkxer('DGEMLQ', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 5;
    dgemlq('L', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W, 1, INFO);
    chkxer('DGEMLQ', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 5;
    dgemlq('R', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W, 1, INFO);
    chkxer('DGEMLQ', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 7;
    dgemlq('L', 'N', 1, 2, 0, A, 0, TAU, 1, C, 1, W, 1, INFO);
    chkxer('DGEMLQ', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 9;
    dgemlq('R', 'N', 2, 2, 1, A, 1, TAU, 0, C, 1, W, 1, INFO);
    chkxer('DGEMLQ', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 9;
    dgemlq('L', 'N', 2, 2, 1, A, 1, TAU, 0, C, 1, W, 1, INFO);
    chkxer('DGEMLQ', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 11;
    dgemlq('L', 'N', 1, 2, 1, A, 1, TAU, 6, C, 0, W, 1, INFO);
    chkxer('DGEMLQ', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 13;
    dgemlq('L', 'N', 2, 2, 1, A, 2, TAU, 6, C, 2, W, 0, INFO);
    chkxer('DGEMLQ', infoc.INFOT, NOUT, LERR, OK, test);
  });

  // Print a summary line.
  alaesm(PATH, OK.value, NOUT);
}
