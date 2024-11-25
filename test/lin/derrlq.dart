// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';

import '../test_driver.dart';
import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void derrlq(final String PATH, final Nout NUNIT, final TestDriver test) {
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

  // Error exits for LQ factorization

  test('DGELQF', () {
    srnamc.SRNAMT = 'DGELQF';
    infoc.INFOT = 1;
    dgelqf(-1, 0, A, 1, B, W, 1, INFO);
    chkxer('DGELQF', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dgelqf(0, -1, A, 1, B, W, 1, INFO);
    chkxer('DGELQF', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 4;
    dgelqf(2, 1, A, 1, B, W, 2, INFO);
    chkxer('DGELQF', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 7;
    dgelqf(2, 1, A, 2, B, W, 1, INFO);
    chkxer('DGELQF', infoc.INFOT, NOUT, LERR, OK);
  });

  test('DGELQ2', () {
    srnamc.SRNAMT = 'DGELQ2';
    infoc.INFOT = 1;
    dgelq2(-1, 0, A, 1, B, W, INFO);
    chkxer('DGELQ2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dgelq2(0, -1, A, 1, B, W, INFO);
    chkxer('DGELQ2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 4;
    dgelq2(2, 1, A, 1, B, W, INFO);
    chkxer('DGELQ2', infoc.INFOT, NOUT, LERR, OK);
  });

  test('DORGLQ', () {
    srnamc.SRNAMT = 'DORGLQ';
    infoc.INFOT = 1;
    dorglq(-1, 0, 0, A, 1, X, W, 1, INFO);
    chkxer('DORGLQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dorglq(0, -1, 0, A, 1, X, W, 1, INFO);
    chkxer('DORGLQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dorglq(2, 1, 0, A, 2, X, W, 2, INFO);
    chkxer('DORGLQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 3;
    dorglq(0, 0, -1, A, 1, X, W, 1, INFO);
    chkxer('DORGLQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 3;
    dorglq(1, 1, 2, A, 1, X, W, 1, INFO);
    chkxer('DORGLQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dorglq(2, 2, 0, A, 1, X, W, 2, INFO);
    chkxer('DORGLQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 8;
    dorglq(2, 2, 0, A, 2, X, W, 1, INFO);
    chkxer('DORGLQ', infoc.INFOT, NOUT, LERR, OK);
  });

  test('DORGL2', () {
    srnamc.SRNAMT = 'DORGL2';
    infoc.INFOT = 1;
    dorgl2(-1, 0, 0, A, 1, X, W, INFO);
    chkxer('DORGL2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dorgl2(0, -1, 0, A, 1, X, W, INFO);
    chkxer('DORGL2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dorgl2(2, 1, 0, A, 2, X, W, INFO);
    chkxer('DORGL2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 3;
    dorgl2(0, 0, -1, A, 1, X, W, INFO);
    chkxer('DORGL2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 3;
    dorgl2(1, 1, 2, A, 1, X, W, INFO);
    chkxer('DORGL2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dorgl2(2, 2, 0, A, 1, X, W, INFO);
    chkxer('DORGL2', infoc.INFOT, NOUT, LERR, OK);
  });

  test('DORMLQ', () {
    srnamc.SRNAMT = 'DORMLQ';
    infoc.INFOT = 1;
    dormlq('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMLQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dormlq('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMLQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 3;
    dormlq('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMLQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 4;
    dormlq('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMLQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dormlq('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMLQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dormlq('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMLQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dormlq('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMLQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 7;
    dormlq('L', 'N', 2, 0, 2, A, 1, X, AF, 2, W, 1, INFO);
    chkxer('DORMLQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 7;
    dormlq('R', 'N', 0, 2, 2, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMLQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 10;
    dormlq('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO);
    chkxer('DORMLQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 12;
    dormlq('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO);
    chkxer('DORMLQ', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 12;
    dormlq('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO);
    chkxer('DORMLQ', infoc.INFOT, NOUT, LERR, OK);
  });

  test('DORML2', () {
    srnamc.SRNAMT = 'DORML2';
    infoc.INFOT = 1;
    dorml2('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO);
    chkxer('DORML2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dorml2('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO);
    chkxer('DORML2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 3;
    dorml2('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO);
    chkxer('DORML2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 4;
    dorml2('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO);
    chkxer('DORML2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dorml2('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO);
    chkxer('DORML2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dorml2('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO);
    chkxer('DORML2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dorml2('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO);
    chkxer('DORML2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 7;
    dorml2('L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, INFO);
    chkxer('DORML2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 7;
    dorml2('R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, INFO);
    chkxer('DORML2', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 10;
    dorml2('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO);
    chkxer('DORML2', infoc.INFOT, NOUT, LERR, OK);
  });

  // Print a summary line.
  alaesm(PATH, OK.value, NOUT);
}
