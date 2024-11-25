// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';

import '../test_driver.dart';
import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void derrge(final String PATH, final Nout NUNIT, final TestDriver test) {
  const NMAX = 4, LW = 3 * NMAX;
  final IP = Array<int>(NMAX), IW = Array<int>(NMAX);
  final A = Matrix<double>(NMAX, NMAX),
      AF = Matrix<double>(NMAX, NMAX),
      B = Array<double>(NMAX),
      R1 = Array<double>(NMAX),
      R2 = Array<double>(NMAX),
      W = Array<double>(LW),
      X = Array<double>(NMAX);
  final INFO = Box(0);
  final ANRM = Box(0.0), CCOND = Box(0.0), RCOND = Box(0.0);

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
    R1[J] = 0.0;
    R2[J] = 0.0;
    W[J] = 0.0;
    X[J] = 0.0;
    IP[J] = J;
    IW[J] = J;
  }
  infoc.OK.value = true;

  if (lsamen(2, C2, 'GE')) {
    // Test error exits of the routines that use the LU decomposition
    // of a general matrix.

    test('DGETRF', () {
      srnamc.SRNAMT = 'DGETRF';
      infoc.INFOT = 1;
      dgetrf(-1, 0, A, 1, IP, INFO);
      chkxer('DGETRF', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgetrf(0, -1, A, 1, IP, INFO);
      chkxer('DGETRF', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dgetrf(2, 1, A, 1, IP, INFO);
      chkxer('DGETRF', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DGETF2', () {
      srnamc.SRNAMT = 'DGETF2';
      infoc.INFOT = 1;
      dgetf2(-1, 0, A, 1, IP, INFO);
      chkxer('DGETF2', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgetf2(0, -1, A, 1, IP, INFO);
      chkxer('DGETF2', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dgetf2(2, 1, A, 1, IP, INFO);
      chkxer('DGETF2', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DGETRI', () {
      srnamc.SRNAMT = 'DGETRI';
      infoc.INFOT = 1;
      dgetri(-1, A, 1, IP, W, LW, INFO);
      chkxer('DGETRI', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dgetri(2, A, 1, IP, W, LW, INFO);
      chkxer('DGETRI', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DGETRS', () {
      srnamc.SRNAMT = 'DGETRS';
      infoc.INFOT = 1;
      dgetrs('/', 0, 0, A, 1, IP, B.asMatrix(), 1, INFO);
      chkxer('DGETRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgetrs('N', -1, 0, A, 1, IP, B.asMatrix(), 1, INFO);
      chkxer('DGETRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dgetrs('N', 0, -1, A, 1, IP, B.asMatrix(), 1, INFO);
      chkxer('DGETRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dgetrs('N', 2, 1, A, 1, IP, B.asMatrix(), 2, INFO);
      chkxer('DGETRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      dgetrs('N', 2, 1, A, 2, IP, B.asMatrix(), 1, INFO);
      chkxer('DGETRS', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DGERFS', () {
      srnamc.SRNAMT = 'DGERFS';
      infoc.INFOT = 1;
      dgerfs('/', 0, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1, R1,
          R2, W, IW, INFO);
      chkxer('DGERFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgerfs('N', -1, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1, R1,
          R2, W, IW, INFO);
      chkxer('DGERFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dgerfs('N', 0, -1, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1, R1,
          R2, W, IW, INFO);
      chkxer('DGERFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dgerfs('N', 2, 1, A, 1, AF, 2, IP, B.asMatrix(), 2, X.asMatrix(), 2, R1,
          R2, W, IW, INFO);
      chkxer('DGERFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 7;
      dgerfs('N', 2, 1, A, 2, AF, 1, IP, B.asMatrix(), 2, X.asMatrix(), 2, R1,
          R2, W, IW, INFO);
      chkxer('DGERFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 10;
      dgerfs('N', 2, 1, A, 2, AF, 2, IP, B.asMatrix(), 1, X.asMatrix(), 2, R1,
          R2, W, IW, INFO);
      chkxer('DGERFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 12;
      dgerfs('N', 2, 1, A, 2, AF, 2, IP, B.asMatrix(), 2, X.asMatrix(), 1, R1,
          R2, W, IW, INFO);
      chkxer('DGERFS', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DGECON', () {
      srnamc.SRNAMT = 'DGECON';
      infoc.INFOT = 1;
      dgecon('/', 0, A, 1, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DGECON', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgecon('1', -1, A, 1, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DGECON', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dgecon('1', 2, A, 1, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DGECON', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DGEEQU', () {
      srnamc.SRNAMT = 'DGEEQU';
      infoc.INFOT = 1;
      dgeequ(-1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
      chkxer('DGEEQU', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgeequ(0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
      chkxer('DGEEQU', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dgeequ(2, 2, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
      chkxer('DGEEQU', infoc.INFOT, NOUT, LERR, OK, test);
    });
  } else if (lsamen(2, C2, 'GB')) {
    // Test error exits of the routines that use the LU decomposition
    // of a general band matrix.

    test('DGBTRF', () {
      srnamc.SRNAMT = 'DGBTRF';
      infoc.INFOT = 1;
      dgbtrf(-1, 0, 0, 0, A, 1, IP, INFO);
      chkxer('DGBTRF', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgbtrf(0, -1, 0, 0, A, 1, IP, INFO);
      chkxer('DGBTRF', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dgbtrf(1, 1, -1, 0, A, 1, IP, INFO);
      chkxer('DGBTRF', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dgbtrf(1, 1, 0, -1, A, 1, IP, INFO);
      chkxer('DGBTRF', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dgbtrf(2, 2, 1, 1, A, 3, IP, INFO);
      chkxer('DGBTRF', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DGBTF2', () {
      srnamc.SRNAMT = 'DGBTF2';
      infoc.INFOT = 1;
      dgbtf2(-1, 0, 0, 0, A, 1, IP, INFO);
      chkxer('DGBTF2', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgbtf2(0, -1, 0, 0, A, 1, IP, INFO);
      chkxer('DGBTF2', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dgbtf2(1, 1, -1, 0, A, 1, IP, INFO);
      chkxer('DGBTF2', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dgbtf2(1, 1, 0, -1, A, 1, IP, INFO);
      chkxer('DGBTF2', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dgbtf2(2, 2, 1, 1, A, 3, IP, INFO);
      chkxer('DGBTF2', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DGBTRS', () {
      srnamc.SRNAMT = 'DGBTRS';
      infoc.INFOT = 1;
      dgbtrs('/', 0, 0, 0, 1, A, 1, IP, B.asMatrix(), 1, INFO);
      chkxer('DGBTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgbtrs('N', -1, 0, 0, 1, A, 1, IP, B.asMatrix(), 1, INFO);
      chkxer('DGBTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dgbtrs('N', 1, -1, 0, 1, A, 1, IP, B.asMatrix(), 1, INFO);
      chkxer('DGBTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dgbtrs('N', 1, 0, -1, 1, A, 1, IP, B.asMatrix(), 1, INFO);
      chkxer('DGBTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dgbtrs('N', 1, 0, 0, -1, A, 1, IP, B.asMatrix(), 1, INFO);
      chkxer('DGBTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 7;
      dgbtrs('N', 2, 1, 1, 1, A, 3, IP, B.asMatrix(), 2, INFO);
      chkxer('DGBTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 10;
      dgbtrs('N', 2, 0, 0, 1, A, 1, IP, B.asMatrix(), 1, INFO);
      chkxer('DGBTRS', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DGBRFS', () {
      srnamc.SRNAMT = 'DGBRFS';
      infoc.INFOT = 1;
      dgbrfs('/', 0, 0, 0, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
          R1, R2, W, IW, INFO);
      chkxer('DGBRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgbrfs('N', -1, 0, 0, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(),
          1, R1, R2, W, IW, INFO);
      chkxer('DGBRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dgbrfs('N', 1, -1, 0, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(),
          1, R1, R2, W, IW, INFO);
      chkxer('DGBRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dgbrfs('N', 1, 0, -1, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(),
          1, R1, R2, W, IW, INFO);
      chkxer('DGBRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dgbrfs('N', 1, 0, 0, -1, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(),
          1, R1, R2, W, IW, INFO);
      chkxer('DGBRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 7;
      dgbrfs('N', 2, 1, 1, 1, A, 2, AF, 4, IP, B.asMatrix(), 2, X.asMatrix(), 2,
          R1, R2, W, IW, INFO);
      chkxer('DGBRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 9;
      dgbrfs('N', 2, 1, 1, 1, A, 3, AF, 3, IP, B.asMatrix(), 2, X.asMatrix(), 2,
          R1, R2, W, IW, INFO);
      chkxer('DGBRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 12;
      dgbrfs('N', 2, 0, 0, 1, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 2,
          R1, R2, W, IW, INFO);
      chkxer('DGBRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 14;
      dgbrfs('N', 2, 0, 0, 1, A, 1, AF, 1, IP, B.asMatrix(), 2, X.asMatrix(), 1,
          R1, R2, W, IW, INFO);
      chkxer('DGBRFS', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DGBCON', () {
      srnamc.SRNAMT = 'DGBCON';
      infoc.INFOT = 1;
      dgbcon('/', 0, 0, 0, A, 1, IP, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DGBCON', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgbcon('1', -1, 0, 0, A, 1, IP, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DGBCON', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dgbcon('1', 1, -1, 0, A, 1, IP, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DGBCON', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dgbcon('1', 1, 0, -1, A, 1, IP, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DGBCON', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dgbcon('1', 2, 1, 1, A, 3, IP, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DGBCON', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DGBEQU', () {
      srnamc.SRNAMT = 'DGBEQU';
      infoc.INFOT = 1;
      dgbequ(-1, 0, 0, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
      chkxer('DGBEQU', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgbequ(0, -1, 0, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
      chkxer('DGBEQU', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dgbequ(1, 1, -1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
      chkxer('DGBEQU', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dgbequ(1, 1, 0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
      chkxer('DGBEQU', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dgbequ(2, 2, 1, 1, A, 2, R1, R2, RCOND, CCOND, ANRM, INFO);
      chkxer('DGBEQU', infoc.INFOT, NOUT, LERR, OK, test);
    });
  }

  // Print a summary line.
  alaesm(PATH, infoc.OK.value, infoc.NOUT);
}
