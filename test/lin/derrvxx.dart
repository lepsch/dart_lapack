// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';

import '../test_driver.dart';
import 'chkxer.dart';
import 'common.dart';

void derrvx(final String PATH, final Nout NUNIT, [final TestDriver? test]) {
  const NMAX = 4;
  const ONE = 1.0;
  int N_ERR_BNDS, NPARAMS;
  final IP = Array<int>(NMAX), IW = Array<int>(NMAX);
  final A = Matrix<double>(NMAX, NMAX),
      AF = Matrix<double>(NMAX, NMAX),
      B = Array<double>(NMAX),
      C = Array<double>(NMAX),
      E = Array<double>(NMAX),
      R = Array<double>(NMAX),
      R1 = Array<double>(NMAX),
      R2 = Array<double>(NMAX),
      W = Array<double>(2 * NMAX),
      X = Array<double>(NMAX),
      ERR_BNDS_N = Matrix<double>(NMAX, 3),
      ERR_BNDS_C = Matrix<double>(NMAX, 3),
      PARAMS = Array<double>(1),
      BERR = Array<double>(1);
  final INFO = Box(0);
  final RCOND = Box(0.0), RPVGRW = Box(0.0);
  final EQ = Box('');

  final NOUT = infoc.NOUT = NUNIT;
  final OK = infoc.OK;
  final LERR = infoc.LERR;
  NOUT.println();
  final C2 = PATH.substring(1, 3);

  // Set the variables to innocuous values.

  for (var J = 1; J <= NMAX; J++) {
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = 1.0 / (I + J);
      AF[I][J] = 1.0 / (I + J);
    }
    B[J] = 0.0;
    E[J] = 0.0;
    R1[J] = 0.0;
    R2[J] = 0.0;
    W[J] = 0.0;
    X[J] = 0.0;
    C[J] = 0.0;
    R[J] = 0.0;
    IP[J] = J;
  }
  EQ.value = ' ';
  OK.value = true;

  if (lsamen(2, C2, 'GE')) {
    test?.test('DGESV', () {
      srnamc.SRNAMT = 'DGESV';
      infoc.INFOT = 1;
      dgesv(-1, 0, A, 1, IP, B.asMatrix(), 1, INFO);
      chkxer('DGESV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgesv(0, -1, A, 1, IP, B.asMatrix(), 1, INFO);
      chkxer('DGESV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dgesv(2, 1, A, 1, IP, B.asMatrix(), 2, INFO);
      chkxer('DGESV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 7;
      dgesv(2, 1, A, 2, IP, B.asMatrix(), 1, INFO);
      chkxer('DGESV', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test?.test('DGESVX', () {
      srnamc.SRNAMT = 'DGESVX';
      infoc.INFOT = 1;
      dgesvx('/', 'N', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGESVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgesvx('N', '/', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGESVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dgesvx('N', 'N', -1, 0, A, 1, AF, 1, IP, EQ, R, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGESVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dgesvx('N', 'N', 0, -1, A, 1, AF, 1, IP, EQ, R, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGESVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dgesvx('N', 'N', 2, 1, A, 1, AF, 2, IP, EQ, R, C, B.asMatrix(), 2,
          X.asMatrix(), 2, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGESVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      dgesvx('N', 'N', 2, 1, A, 2, AF, 1, IP, EQ, R, C, B.asMatrix(), 2,
          X.asMatrix(), 2, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGESVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 10;
      EQ.value = '/';
      dgesvx('F', 'N', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGESVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 11;
      EQ.value = 'R';
      dgesvx('F', 'N', 1, 0, A, 1, AF, 1, IP, EQ, R, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGESVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 12;
      EQ.value = 'C';
      dgesvx('F', 'N', 1, 0, A, 1, AF, 1, IP, EQ, R, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGESVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 14;
      dgesvx('N', 'N', 2, 1, A, 2, AF, 2, IP, EQ, R, C, B.asMatrix(), 1,
          X.asMatrix(), 2, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGESVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 16;
      dgesvx('N', 'N', 2, 1, A, 2, AF, 2, IP, EQ, R, C, B.asMatrix(), 2,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGESVX', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test?.test('DGESVXX', () {
      N_ERR_BNDS = 3;
      NPARAMS = 1;
      srnamc.SRNAMT = 'DGESVXX';
      infoc.INFOT = 1;
      dgesvxx(
          '/',
          'N',
          0,
          0,
          A,
          1,
          AF,
          1,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGESVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgesvxx(
          'N',
          '/',
          0,
          0,
          A,
          1,
          AF,
          1,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGESVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dgesvxx(
          'N',
          'N',
          -1,
          0,
          A,
          1,
          AF,
          1,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGESVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dgesvxx(
          'N',
          'N',
          0,
          -1,
          A,
          1,
          AF,
          1,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGESVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dgesvxx(
          'N',
          'N',
          2,
          1,
          A,
          1,
          AF,
          2,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          2,
          X.asMatrix(),
          2,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGESVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      dgesvxx(
          'N',
          'N',
          2,
          1,
          A,
          2,
          AF,
          1,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          2,
          X.asMatrix(),
          2,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGESVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 10;
      EQ.value = '/';
      dgesvxx(
          'F',
          'N',
          0,
          0,
          A,
          1,
          AF,
          1,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGESVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 11;
      EQ.value = 'R';
      dgesvxx(
          'F',
          'N',
          1,
          0,
          A,
          1,
          AF,
          1,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGESVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 12;
      EQ.value = 'C';
      dgesvxx(
          'F',
          'N',
          1,
          0,
          A,
          1,
          AF,
          1,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGESVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 14;
      dgesvxx(
          'N',
          'N',
          2,
          1,
          A,
          2,
          AF,
          2,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          2,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGESVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 16;
      dgesvxx(
          'N',
          'N',
          2,
          1,
          A,
          2,
          AF,
          2,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          2,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGESVXX', infoc.INFOT, NOUT, LERR, OK, test);
    });
  } else if (lsamen(2, C2, 'GB')) {
    test?.test('DGBSV', () {
      srnamc.SRNAMT = 'DGBSV';
      infoc.INFOT = 1;
      dgbsv(-1, 0, 0, 0, A, 1, IP, B.asMatrix(), 1, INFO);
      chkxer('DGBSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgbsv(1, -1, 0, 0, A, 1, IP, B.asMatrix(), 1, INFO);
      chkxer('DGBSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dgbsv(1, 0, -1, 0, A, 1, IP, B.asMatrix(), 1, INFO);
      chkxer('DGBSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dgbsv(0, 0, 0, -1, A, 1, IP, B.asMatrix(), 1, INFO);
      chkxer('DGBSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dgbsv(1, 1, 1, 0, A, 3, IP, B.asMatrix(), 1, INFO);
      chkxer('DGBSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 9;
      dgbsv(2, 0, 0, 0, A, 1, IP, B.asMatrix(), 1, INFO);
      chkxer('DGBSV', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test?.test('DGBSVX', () {
      srnamc.SRNAMT = 'DGBSVX';
      infoc.INFOT = 1;
      dgbsvx('/', 'N', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGBSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgbsvx('N', '/', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGBSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dgbsvx('N', 'N', -1, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGBSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dgbsvx('N', 'N', 1, -1, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGBSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dgbsvx('N', 'N', 1, 0, -1, 0, A, 1, AF, 1, IP, EQ, R, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGBSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dgbsvx('N', 'N', 0, 0, 0, -1, A, 1, AF, 1, IP, EQ, R, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGBSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      dgbsvx('N', 'N', 1, 1, 1, 0, A, 2, AF, 4, IP, EQ, R, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGBSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 10;
      dgbsvx('N', 'N', 1, 1, 1, 0, A, 3, AF, 3, IP, EQ, R, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGBSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 12;
      EQ.value = '/';
      dgbsvx('F', 'N', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGBSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 13;
      EQ.value = 'R';
      dgbsvx('F', 'N', 1, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGBSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 14;
      EQ.value = 'C';
      dgbsvx('F', 'N', 1, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGBSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 16;
      dgbsvx('N', 'N', 2, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B.asMatrix(), 1,
          X.asMatrix(), 2, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGBSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 18;
      dgbsvx('N', 'N', 2, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B.asMatrix(), 2,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DGBSVX', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test?.test('DGBSVXX', () {
      N_ERR_BNDS = 3;
      NPARAMS = 1;
      srnamc.SRNAMT = 'DGBSVXX';
      infoc.INFOT = 1;
      dgbsvxx(
          '/',
          'N',
          0,
          0,
          0,
          0,
          A,
          1,
          AF,
          1,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGBSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgbsvxx(
          'N',
          '/',
          0,
          1,
          1,
          0,
          A,
          1,
          AF,
          1,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGBSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dgbsvxx(
          'N',
          'N',
          -1,
          1,
          1,
          0,
          A,
          1,
          AF,
          1,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGBSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dgbsvxx(
          'N',
          'N',
          2,
          -1,
          1,
          0,
          A,
          1,
          AF,
          1,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGBSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dgbsvxx(
          'N',
          'N',
          2,
          1,
          -1,
          0,
          A,
          1,
          AF,
          1,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGBSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dgbsvxx(
          'N',
          'N',
          0,
          1,
          1,
          -1,
          A,
          1,
          AF,
          1,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGBSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      dgbsvxx(
          'N',
          'N',
          2,
          1,
          1,
          1,
          A,
          2,
          AF,
          2,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          2,
          X.asMatrix(),
          2,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGBSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 10;
      dgbsvxx(
          'N',
          'N',
          2,
          1,
          1,
          1,
          A,
          3,
          AF,
          3,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          2,
          X.asMatrix(),
          2,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGBSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 12;
      EQ.value = '/';
      dgbsvxx(
          'F',
          'N',
          0,
          1,
          1,
          0,
          A,
          3,
          AF,
          4,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGBSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 13;
      EQ.value = 'R';
      dgbsvxx(
          'F',
          'N',
          1,
          1,
          1,
          0,
          A,
          3,
          AF,
          4,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGBSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 14;
      EQ.value = 'C';
      dgbsvxx(
          'F',
          'N',
          1,
          1,
          1,
          0,
          A,
          3,
          AF,
          4,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGBSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 15;
      dgbsvxx(
          'N',
          'N',
          2,
          1,
          1,
          1,
          A,
          3,
          AF,
          4,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          2,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGBSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 16;
      dgbsvxx(
          'N',
          'N',
          2,
          1,
          1,
          1,
          A,
          3,
          AF,
          4,
          IP,
          EQ,
          R,
          C,
          B.asMatrix(),
          2,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DGBSVXX', infoc.INFOT, NOUT, LERR, OK, test);
    });
  } else if (lsamen(2, C2, 'GT')) {
    test?.test('DGTSV', () {
      srnamc.SRNAMT = 'DGTSV';
      infoc.INFOT = 1;
      dgtsv(-1, 0, A(1, 1).asArray(), A(1, 2).asArray(), A(1, 3).asArray(),
          B.asMatrix(), 1, INFO);
      chkxer('DGTSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgtsv(0, -1, A(1, 1).asArray(), A(1, 2).asArray(), A(1, 3).asArray(),
          B.asMatrix(), 1, INFO);
      chkxer('DGTSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 7;
      dgtsv(2, 0, A(1, 1).asArray(), A(1, 2).asArray(), A(1, 3).asArray(),
          B.asMatrix(), 1, INFO);
      chkxer('DGTSV', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test?.test('DGTSVX', () {
      srnamc.SRNAMT = 'DGTSVX';
      infoc.INFOT = 1;
      dgtsvx(
          '/',
          'N',
          0,
          0,
          A(1, 1).asArray(),
          A(1, 2).asArray(),
          A(1, 3).asArray(),
          AF(1, 1).asArray(),
          AF(1, 2).asArray(),
          AF(1, 3).asArray(),
          AF(1, 4).asArray(),
          IP,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          R1,
          R2,
          W,
          IW,
          INFO);
      chkxer('DGTSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgtsvx(
          'N',
          '/',
          0,
          0,
          A(1, 1).asArray(),
          A(1, 2).asArray(),
          A(1, 3).asArray(),
          AF(1, 1).asArray(),
          AF(1, 2).asArray(),
          AF(1, 3).asArray(),
          AF(1, 4).asArray(),
          IP,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          R1,
          R2,
          W,
          IW,
          INFO);
      chkxer('DGTSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dgtsvx(
          'N',
          'N',
          -1,
          0,
          A(1, 1).asArray(),
          A(1, 2).asArray(),
          A(1, 3).asArray(),
          AF(1, 1).asArray(),
          AF(1, 2).asArray(),
          AF(1, 3).asArray(),
          AF(1, 4).asArray(),
          IP,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          R1,
          R2,
          W,
          IW,
          INFO);
      chkxer('DGTSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dgtsvx(
          'N',
          'N',
          0,
          -1,
          A(1, 1).asArray(),
          A(1, 2).asArray(),
          A(1, 3).asArray(),
          AF(1, 1).asArray(),
          AF(1, 2).asArray(),
          AF(1, 3).asArray(),
          AF(1, 4).asArray(),
          IP,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          R1,
          R2,
          W,
          IW,
          INFO);
      chkxer('DGTSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 14;
      dgtsvx(
          'N',
          'N',
          2,
          0,
          A(1, 1).asArray(),
          A(1, 2).asArray(),
          A(1, 3).asArray(),
          AF(1, 1).asArray(),
          AF(1, 2).asArray(),
          AF(1, 3).asArray(),
          AF(1, 4).asArray(),
          IP,
          B.asMatrix(),
          1,
          X.asMatrix(),
          2,
          RCOND,
          R1,
          R2,
          W,
          IW,
          INFO);
      chkxer('DGTSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 16;
      dgtsvx(
          'N',
          'N',
          2,
          0,
          A(1, 1).asArray(),
          A(1, 2).asArray(),
          A(1, 3).asArray(),
          AF(1, 1).asArray(),
          AF(1, 2).asArray(),
          AF(1, 3).asArray(),
          AF(1, 4).asArray(),
          IP,
          B.asMatrix(),
          2,
          X.asMatrix(),
          1,
          RCOND,
          R1,
          R2,
          W,
          IW,
          INFO);
      chkxer('DGTSVX', infoc.INFOT, NOUT, LERR, OK, test);
    });
  } else if (lsamen(2, C2, 'PO')) {
    test?.test('DPOSV', () {
      srnamc.SRNAMT = 'DPOSV';
      infoc.INFOT = 1;
      dposv('/', 0, 0, A, 1, B.asMatrix(), 1, INFO);
      chkxer('DPOSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dposv('U', -1, 0, A, 1, B.asMatrix(), 1, INFO);
      chkxer('DPOSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dposv('U', 0, -1, A, 1, B.asMatrix(), 1, INFO);
      chkxer('DPOSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dposv('U', 2, 0, A, 1, B.asMatrix(), 2, INFO);
      chkxer('DPOSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 7;
      dposv('U', 2, 0, A, 2, B.asMatrix(), 1, INFO);
      chkxer('DPOSV', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test?.test('DPOSVX', () {
      srnamc.SRNAMT = 'DPOSVX';
      infoc.INFOT = 1;
      dposvx('/', 'U', 0, 0, A, 1, AF, 1, EQ, C, B.asMatrix(), 1, X.asMatrix(),
          1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPOSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dposvx('N', '/', 0, 0, A, 1, AF, 1, EQ, C, B.asMatrix(), 1, X.asMatrix(),
          1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPOSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dposvx('N', 'U', -1, 0, A, 1, AF, 1, EQ, C, B.asMatrix(), 1, X.asMatrix(),
          1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPOSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dposvx('N', 'U', 0, -1, A, 1, AF, 1, EQ, C, B.asMatrix(), 1, X.asMatrix(),
          1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPOSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dposvx('N', 'U', 2, 0, A, 1, AF, 2, EQ, C, B.asMatrix(), 2, X.asMatrix(),
          2, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPOSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      dposvx('N', 'U', 2, 0, A, 2, AF, 1, EQ, C, B.asMatrix(), 2, X.asMatrix(),
          2, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPOSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 9;
      EQ.value = '/';
      dposvx('F', 'U', 0, 0, A, 1, AF, 1, EQ, C, B.asMatrix(), 1, X.asMatrix(),
          1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPOSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 10;
      EQ.value = 'Y';
      dposvx('F', 'U', 1, 0, A, 1, AF, 1, EQ, C, B.asMatrix(), 1, X.asMatrix(),
          1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPOSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 12;
      dposvx('N', 'U', 2, 0, A, 2, AF, 2, EQ, C, B.asMatrix(), 1, X.asMatrix(),
          2, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPOSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 14;
      dposvx('N', 'U', 2, 0, A, 2, AF, 2, EQ, C, B.asMatrix(), 2, X.asMatrix(),
          1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPOSVX', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test?.test('DPOSVXX', () {
      N_ERR_BNDS = 3;
      NPARAMS = 1;
      srnamc.SRNAMT = 'DPOSVXX';
      infoc.INFOT = 1;
      dposvxx(
          '/',
          'U',
          0,
          0,
          A,
          1,
          AF,
          1,
          EQ,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DPOSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dposvxx(
          'N',
          '/',
          0,
          0,
          A,
          1,
          AF,
          1,
          EQ,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DPOSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dposvxx(
          'N',
          'U',
          -1,
          0,
          A,
          1,
          AF,
          1,
          EQ,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DPOSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dposvxx(
          'N',
          'U',
          0,
          -1,
          A,
          1,
          AF,
          1,
          EQ,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DPOSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dposvxx(
          'N',
          'U',
          2,
          0,
          A,
          1,
          AF,
          2,
          EQ,
          C,
          B.asMatrix(),
          2,
          X.asMatrix(),
          2,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DPOSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      dposvxx(
          'N',
          'U',
          2,
          0,
          A,
          2,
          AF,
          1,
          EQ,
          C,
          B.asMatrix(),
          2,
          X.asMatrix(),
          2,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DPOSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 9;
      EQ.value = '/';
      dposvxx(
          'F',
          'U',
          0,
          0,
          A,
          1,
          AF,
          1,
          EQ,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DPOSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 10;
      EQ.value = 'Y';
      dposvxx(
          'F',
          'U',
          1,
          0,
          A,
          1,
          AF,
          1,
          EQ,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DPOSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 12;
      dposvxx(
          'N',
          'U',
          2,
          0,
          A,
          2,
          AF,
          2,
          EQ,
          C,
          B.asMatrix(),
          1,
          X.asMatrix(),
          2,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DPOSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 14;
      dposvxx(
          'N',
          'U',
          2,
          0,
          A,
          2,
          AF,
          2,
          EQ,
          C,
          B.asMatrix(),
          2,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DPOSVXX', infoc.INFOT, NOUT, LERR, OK, test);
    });
  } else if (lsamen(2, C2, 'PP')) {
    test?.test('DPPSV', () {
      srnamc.SRNAMT = 'DPPSV';
      infoc.INFOT = 1;
      dppsv('/', 0, 0, A.asArray(), B.asMatrix(), 1, INFO);
      chkxer('DPPSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dppsv('U', -1, 0, A.asArray(), B.asMatrix(), 1, INFO);
      chkxer('DPPSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dppsv('U', 0, -1, A.asArray(), B.asMatrix(), 1, INFO);
      chkxer('DPPSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dppsv('U', 2, 0, A.asArray(), B.asMatrix(), 1, INFO);
      chkxer('DPPSV', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test?.test('DPPSVX', () {
      srnamc.SRNAMT = 'DPPSVX';
      infoc.INFOT = 1;
      dppsvx('/', 'U', 0, 0, A.asArray(), AF.asArray(), EQ, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPPSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dppsvx('N', '/', 0, 0, A.asArray(), AF.asArray(), EQ, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPPSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dppsvx('N', 'U', -1, 0, A.asArray(), AF.asArray(), EQ, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPPSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dppsvx('N', 'U', 0, -1, A.asArray(), AF.asArray(), EQ, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPPSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 7;
      EQ.value = '/';
      dppsvx('F', 'U', 0, 0, A.asArray(), AF.asArray(), EQ, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPPSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      EQ.value = 'Y';
      dppsvx('F', 'U', 1, 0, A.asArray(), AF.asArray(), EQ, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPPSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 10;
      dppsvx('N', 'U', 2, 0, A.asArray(), AF.asArray(), EQ, C, B.asMatrix(), 1,
          X.asMatrix(), 2, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPPSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 12;
      dppsvx('N', 'U', 2, 0, A.asArray(), AF.asArray(), EQ, C, B.asMatrix(), 2,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPPSVX', infoc.INFOT, NOUT, LERR, OK, test);
    });
  } else if (lsamen(2, C2, 'PB')) {
    test?.test('DPBSV', () {
      srnamc.SRNAMT = 'DPBSV';
      infoc.INFOT = 1;
      dpbsv('/', 0, 0, 0, A, 1, B.asMatrix(), 1, INFO);
      chkxer('DPBSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dpbsv('U', -1, 0, 0, A, 1, B.asMatrix(), 1, INFO);
      chkxer('DPBSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dpbsv('U', 1, -1, 0, A, 1, B.asMatrix(), 1, INFO);
      chkxer('DPBSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dpbsv('U', 0, 0, -1, A, 1, B.asMatrix(), 1, INFO);
      chkxer('DPBSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dpbsv('U', 1, 1, 0, A, 1, B.asMatrix(), 2, INFO);
      chkxer('DPBSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      dpbsv('U', 2, 0, 0, A, 1, B.asMatrix(), 1, INFO);
      chkxer('DPBSV', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test?.test('DPBSVX', () {
      srnamc.SRNAMT = 'DPBSVX';
      infoc.INFOT = 1;
      dpbsvx('/', 'U', 0, 0, 0, A, 1, AF, 1, EQ, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPBSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dpbsvx('N', '/', 0, 0, 0, A, 1, AF, 1, EQ, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPBSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dpbsvx('N', 'U', -1, 0, 0, A, 1, AF, 1, EQ, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPBSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dpbsvx('N', 'U', 1, -1, 0, A, 1, AF, 1, EQ, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPBSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dpbsvx('N', 'U', 0, 0, -1, A, 1, AF, 1, EQ, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPBSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 7;
      dpbsvx('N', 'U', 1, 1, 0, A, 1, AF, 2, EQ, C, B.asMatrix(), 2,
          X.asMatrix(), 2, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPBSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 9;
      dpbsvx('N', 'U', 1, 1, 0, A, 2, AF, 1, EQ, C, B.asMatrix(), 2,
          X.asMatrix(), 2, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPBSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 10;
      EQ.value = '/';
      dpbsvx('F', 'U', 0, 0, 0, A, 1, AF, 1, EQ, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPBSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 11;
      EQ.value = 'Y';
      dpbsvx('F', 'U', 1, 0, 0, A, 1, AF, 1, EQ, C, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPBSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 13;
      dpbsvx('N', 'U', 2, 0, 0, A, 1, AF, 1, EQ, C, B.asMatrix(), 1,
          X.asMatrix(), 2, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPBSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 15;
      dpbsvx('N', 'U', 2, 0, 0, A, 1, AF, 1, EQ, C, B.asMatrix(), 2,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DPBSVX', infoc.INFOT, NOUT, LERR, OK, test);
    });
  } else if (lsamen(2, C2, 'PT')) {
    test?.test('DPTSV', () {
      srnamc.SRNAMT = 'DPTSV';
      infoc.INFOT = 1;
      dptsv(-1, 0, A(1, 1).asArray(), A(1, 2).asArray(), B.asMatrix(), 1, INFO);
      chkxer('DPTSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dptsv(0, -1, A(1, 1).asArray(), A(1, 2).asArray(), B.asMatrix(), 1, INFO);
      chkxer('DPTSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dptsv(2, 0, A(1, 1).asArray(), A(1, 2).asArray(), B.asMatrix(), 1, INFO);
      chkxer('DPTSV', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test?.test('DPTSVX', () {
      srnamc.SRNAMT = 'DPTSVX';
      infoc.INFOT = 1;
      dptsvx(
          '/',
          0,
          0,
          A(1, 1).asArray(),
          A(1, 2).asArray(),
          AF(1, 1).asArray(),
          AF(1, 2).asArray(),
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          R1,
          R2,
          W,
          INFO);
      chkxer('DPTSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dptsvx(
          'N',
          -1,
          0,
          A(1, 1).asArray(),
          A(1, 2).asArray(),
          AF(1, 1).asArray(),
          AF(1, 2).asArray(),
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          R1,
          R2,
          W,
          INFO);
      chkxer('DPTSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dptsvx(
          'N',
          0,
          -1,
          A(1, 1).asArray(),
          A(1, 2).asArray(),
          AF(1, 1).asArray(),
          AF(1, 2).asArray(),
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          R1,
          R2,
          W,
          INFO);
      chkxer('DPTSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 9;
      dptsvx(
          'N',
          2,
          0,
          A(1, 1).asArray(),
          A(1, 2).asArray(),
          AF(1, 1).asArray(),
          AF(1, 2).asArray(),
          B.asMatrix(),
          1,
          X.asMatrix(),
          2,
          RCOND,
          R1,
          R2,
          W,
          INFO);
      chkxer('DPTSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 11;
      dptsvx(
          'N',
          2,
          0,
          A(1, 1).asArray(),
          A(1, 2).asArray(),
          AF(1, 1).asArray(),
          AF(1, 2).asArray(),
          B.asMatrix(),
          2,
          X.asMatrix(),
          1,
          RCOND,
          R1,
          R2,
          W,
          INFO);
      chkxer('DPTSVX', infoc.INFOT, NOUT, LERR, OK, test);
    });
  } else if (lsamen(2, C2, 'SY')) {
    test?.test('DSYSV', () {
      srnamc.SRNAMT = 'DSYSV';
      infoc.INFOT = 1;
      dsysv('/', 0, 0, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
      chkxer('DSYSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dsysv('U', -1, 0, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
      chkxer('DSYSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dsysv('U', 0, -1, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
      chkxer('DSYSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dsysv('U', 2, 0, A, 1, IP, B.asMatrix(), 2, W, 1, INFO);
      chkxer('DSYSV_ROOK', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      dsysv('U', 2, 0, A, 2, IP, B.asMatrix(), 1, W, 1, INFO);
      chkxer('DSYSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 10;
      dsysv('U', 0, 0, A, 1, IP, B.asMatrix(), 1, W, 0, INFO);
      chkxer('DSYSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 10;
      dsysv('U', 0, 0, A, 1, IP, B.asMatrix(), 1, W, -2, INFO);
      chkxer('DSYSV', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test?.test('DSYSVX', () {
      srnamc.SRNAMT = 'DSYSVX';
      infoc.INFOT = 1;
      dsysvx('/', 'U', 0, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
          RCOND, R1, R2, W, 1, IW, INFO);
      chkxer('DSYSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dsysvx('N', '/', 0, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
          RCOND, R1, R2, W, 1, IW, INFO);
      chkxer('DSYSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dsysvx('N', 'U', -1, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
          RCOND, R1, R2, W, 1, IW, INFO);
      chkxer('DSYSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dsysvx('N', 'U', 0, -1, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
          RCOND, R1, R2, W, 1, IW, INFO);
      chkxer('DSYSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dsysvx('N', 'U', 2, 0, A, 1, AF, 2, IP, B.asMatrix(), 2, X.asMatrix(), 2,
          RCOND, R1, R2, W, 4, IW, INFO);
      chkxer('DSYSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      dsysvx('N', 'U', 2, 0, A, 2, AF, 1, IP, B.asMatrix(), 2, X.asMatrix(), 2,
          RCOND, R1, R2, W, 4, IW, INFO);
      chkxer('DSYSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 11;
      dsysvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B.asMatrix(), 1, X.asMatrix(), 2,
          RCOND, R1, R2, W, 4, IW, INFO);
      chkxer('DSYSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 13;
      dsysvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B.asMatrix(), 2, X.asMatrix(), 1,
          RCOND, R1, R2, W, 4, IW, INFO);
      chkxer('DSYSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 18;
      dsysvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B.asMatrix(), 2, X.asMatrix(), 2,
          RCOND, R1, R2, W, 3, IW, INFO);
      chkxer('DSYSVX', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test?.test('DSYSVXX', () {
      N_ERR_BNDS = 3;
      NPARAMS = 1;
      srnamc.SRNAMT = 'DSYSVXX';
      infoc.INFOT = 1;
      EQ.value = 'N';
      dsysvxx(
          '/',
          'U',
          0,
          0,
          A,
          1,
          AF,
          1,
          IP,
          EQ,
          R,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DSYSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dsysvxx(
          'N',
          '/',
          0,
          0,
          A,
          1,
          AF,
          1,
          IP,
          EQ,
          R,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DSYSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dsysvxx(
          'N',
          'U',
          -1,
          0,
          A,
          1,
          AF,
          1,
          IP,
          EQ,
          R,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DSYSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      EQ.value = '/';
      dsysvxx(
          'N',
          'U',
          0,
          -1,
          A,
          1,
          AF,
          1,
          IP,
          EQ,
          R,
          B.asMatrix(),
          1,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DSYSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      EQ.value = 'Y';
      infoc.INFOT = 6;
      dsysvxx(
          'N',
          'U',
          2,
          0,
          A,
          1,
          AF,
          2,
          IP,
          EQ,
          R,
          B.asMatrix(),
          2,
          X.asMatrix(),
          2,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DSYSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      dsysvxx(
          'N',
          'U',
          2,
          0,
          A,
          2,
          AF,
          1,
          IP,
          EQ,
          R,
          B.asMatrix(),
          2,
          X.asMatrix(),
          2,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DSYSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 10;
      dsysvxx(
          'F',
          'U',
          2,
          0,
          A,
          2,
          AF,
          2,
          IP,
          Box('A'),
          R,
          B.asMatrix(),
          2,
          X.asMatrix(),
          2,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DSYSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 11;
      EQ.value = 'Y';
      dsysvxx(
          'F',
          'U',
          2,
          0,
          A,
          2,
          AF,
          2,
          IP,
          EQ,
          R,
          B.asMatrix(),
          2,
          X.asMatrix(),
          2,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DSYSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 11;
      EQ.value = 'Y';
      R[1] = -ONE;
      dsysvxx(
          'F',
          'U',
          2,
          0,
          A,
          2,
          AF,
          2,
          IP,
          EQ,
          R,
          B.asMatrix(),
          2,
          X.asMatrix(),
          2,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DSYSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 13;
      EQ.value = 'N';
      dsysvxx(
          'N',
          'U',
          2,
          0,
          A,
          2,
          AF,
          2,
          IP,
          EQ,
          R,
          B.asMatrix(),
          1,
          X.asMatrix(),
          2,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DSYSVXX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 15;
      dsysvxx(
          'N',
          'U',
          2,
          0,
          A,
          2,
          AF,
          2,
          IP,
          EQ,
          R,
          B.asMatrix(),
          2,
          X.asMatrix(),
          1,
          RCOND,
          RPVGRW,
          BERR,
          N_ERR_BNDS,
          ERR_BNDS_N,
          ERR_BNDS_C,
          NPARAMS,
          PARAMS,
          W,
          IW,
          INFO);
      chkxer('DSYSVXX', infoc.INFOT, NOUT, LERR, OK, test);
    });
  } else if (lsamen(2, C2, 'SR')) {
    test?.test('DSYSV_ROOK', () {
      srnamc.SRNAMT = 'DSYSV_ROOK';
      infoc.INFOT = 1;
      dsysv_rook('/', 0, 0, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
      chkxer('DSYSV_ROOK', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dsysv_rook('U', -1, 0, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
      chkxer('DSYSV_ROOK', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dsysv_rook('U', 0, -1, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
      chkxer('DSYSV_ROOK', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dsysv_rook('U', 2, 0, A, 1, IP, B.asMatrix(), 2, W, 1, INFO);
      chkxer('DSYSV_ROOK', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      dsysv_rook('U', 2, 0, A, 2, IP, B.asMatrix(), 1, W, 1, INFO);
      chkxer('DSYSV_ROOK', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 10;
      dsysv_rook('U', 0, 0, A, 1, IP, B.asMatrix(), 1, W, 0, INFO);
      chkxer('DSYSV_ROOK', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 10;
      dsysv_rook('U', 0, 0, A, 1, IP, B.asMatrix(), 1, W, -2, INFO);
      chkxer('DSYSV_ROOK', infoc.INFOT, NOUT, LERR, OK, test);
    });
  } else if (lsamen(2, C2, 'SK')) {
    test?.test('DSYSV_RK', () {
      // Test error exits of the driver that uses factorization
      // of a symmetric indefinite matrix with rook
      // (bounded Bunch-Kaufman) pivoting with the new storage
      // format for factors L ( or U) and D.

      // L (or U) is stored in A, diagonal of D is stored on the
      // diagonal of A, subdiagonal of D is stored in a separate array E.

      srnamc.SRNAMT = 'DSYSV_RK';
      infoc.INFOT = 1;
      dsysv_rk('/', 0, 0, A, 1, E, IP, B.asMatrix(), 1, W, 1, INFO);
      chkxer('DSYSV_RK', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dsysv_rk('U', -1, 0, A, 1, E, IP, B.asMatrix(), 1, W, 1, INFO);
      chkxer('DSYSV_RK', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dsysv_rk('U', 0, -1, A, 1, E, IP, B.asMatrix(), 1, W, 1, INFO);
      chkxer('DSYSV_RK', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dsysv_rk('U', 2, 0, A, 1, E, IP, B.asMatrix(), 2, W, 1, INFO);
      chkxer('DSYSV_RK', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 9;
      dsysv_rk('U', 2, 0, A, 2, E, IP, B.asMatrix(), 1, W, 1, INFO);
      chkxer('DSYSV_RK', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 11;
      dsysv_rk('U', 0, 0, A, 1, E, IP, B.asMatrix(), 1, W, 0, INFO);
      chkxer('DSYSV_RK', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 11;
      dsysv_rk('U', 0, 0, A, 1, E, IP, B.asMatrix(), 1, W, -2, INFO);
      chkxer('DSYSV_RK', infoc.INFOT, NOUT, LERR, OK, test);
    });
  } else if (lsamen(2, C2, 'SP')) {
    test?.test('DSPSV', () {
      srnamc.SRNAMT = 'DSPSV';
      infoc.INFOT = 1;
      dspsv('/', 0, 0, A.asArray(), IP, B.asMatrix(), 1, INFO);
      chkxer('DSPSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dspsv('U', -1, 0, A.asArray(), IP, B.asMatrix(), 1, INFO);
      chkxer('DSPSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dspsv('U', 0, -1, A.asArray(), IP, B.asMatrix(), 1, INFO);
      chkxer('DSPSV', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 7;
      dspsv('U', 2, 0, A.asArray(), IP, B.asMatrix(), 1, INFO);
      chkxer('DSPSV', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test?.test('DSPSVX', () {
      srnamc.SRNAMT = 'DSPSVX';
      infoc.INFOT = 1;
      dspsvx('/', 'U', 0, 0, A.asArray(), AF.asArray(), IP, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DSPSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dspsvx('N', '/', 0, 0, A.asArray(), AF.asArray(), IP, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DSPSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dspsvx('N', 'U', -1, 0, A.asArray(), AF.asArray(), IP, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DSPSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dspsvx('N', 'U', 0, -1, A.asArray(), AF.asArray(), IP, B.asMatrix(), 1,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DSPSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 9;
      dspsvx('N', 'U', 2, 0, A.asArray(), AF.asArray(), IP, B.asMatrix(), 1,
          X.asMatrix(), 2, RCOND, R1, R2, W, IW, INFO);
      chkxer('DSPSVX', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 11;
      dspsvx('N', 'U', 2, 0, A.asArray(), AF.asArray(), IP, B.asMatrix(), 2,
          X.asMatrix(), 1, RCOND, R1, R2, W, IW, INFO);
      chkxer('DSPSVX', infoc.INFOT, NOUT, LERR, OK, test);
    });
  }

  // Print a summary line.

  if (infoc.OK.value) {
    infoc.NOUT
        .println(' ${PATH.a3} drivers passed the tests of the error exits');
  } else {
    NOUT.println(
        ' *** ${PATH.a3} drivers failed the tests of the error exits ***');
  }
}
