// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/ztplqt.dart';
import 'package:dart_lapack/src/ztplqt2.dart';
import 'package:dart_lapack/src/ztpmlqt.dart';

import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void zerrlqtp(final String PATH, final Nout NUNIT) {
  const NMAX = 2;
  final A = Matrix<Complex>(NMAX, NMAX),
      T = Matrix<Complex>(NMAX, NMAX),
      W = Array<Complex>(NMAX),
      B = Matrix<Complex>(NMAX, NMAX),
      C = Matrix<Complex>(NMAX, NMAX);
  final INFO = Box(0);

  final NOUT = NUNIT;
  NOUT.println();

  // Set the variables to innocuous values.

  for (var J = 1; J <= NMAX; J++) {
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = Complex.one / (I + J).toComplex();
      C[I][J] = Complex.one / (I + J).toComplex();
      T[I][J] = Complex.one / (I + J).toComplex();
    }
    W[J] = Complex.zero;
  }
  infoc.OK.value = true;

  // Error exits for TPLQT factorization

  // ZTPLQT

  srnamc.SRNAMT = 'ZTPLQT';
  infoc.INFOT = 1;
  ztplqt(-1, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO);
  chkxer('ZTPLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  ztplqt(1, -1, 0, 1, A, 1, B, 1, T, 1, W, INFO);
  chkxer('ZTPLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  ztplqt(0, 1, -1, 1, A, 1, B, 1, T, 1, W, INFO);
  chkxer('ZTPLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  ztplqt(0, 1, 1, 1, A, 1, B, 1, T, 1, W, INFO);
  chkxer('ZTPLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  ztplqt(0, 1, 0, 0, A, 1, B, 1, T, 1, W, INFO);
  chkxer('ZTPLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  ztplqt(1, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO);
  chkxer('ZTPLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 6;
  ztplqt(2, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO);
  chkxer('ZTPLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  ztplqt(2, 1, 0, 1, A, 2, B, 1, T, 1, W, INFO);
  chkxer('ZTPLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 10;
  ztplqt(2, 2, 1, 2, A, 2, B, 2, T, 1, W, INFO);
  chkxer('ZTPLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // ZTPLQT2

  srnamc.SRNAMT = 'ZTPLQT2';
  infoc.INFOT = 1;
  ztplqt2(-1, 0, 0, A, 1, B, 1, T, 1, INFO);
  chkxer('ZTPLQT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  ztplqt2(0, -1, 0, A, 1, B, 1, T, 1, INFO);
  chkxer('ZTPLQT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  ztplqt2(0, 0, -1, A, 1, B, 1, T, 1, INFO);
  chkxer('ZTPLQT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  ztplqt2(2, 2, 0, A, 1, B, 2, T, 2, INFO);
  chkxer('ZTPLQT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  ztplqt2(2, 2, 0, A, 2, B, 1, T, 2, INFO);
  chkxer('ZTPLQT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 9;
  ztplqt2(2, 2, 0, A, 2, B, 2, T, 1, INFO);
  chkxer('ZTPLQT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // ZTPMLQT

  srnamc.SRNAMT = 'ZTPMLQT';
  infoc.INFOT = 1;
  ztpmlqt('/', 'N', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
  chkxer('ZTPMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  ztpmlqt('L', '/', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
  chkxer('ZTPMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  ztpmlqt('L', 'N', -1, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
  chkxer('ZTPMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  ztpmlqt('L', 'N', 0, -1, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
  chkxer('ZTPMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  ztpmlqt('L', 'N', 0, 0, -1, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
  infoc.INFOT = 6;
  ztpmlqt('L', 'N', 0, 0, 0, -1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
  chkxer('ZTPMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  ztpmlqt('L', 'N', 0, 0, 0, 0, 0, A, 1, T, 1, B, 1, C, 1, W, INFO);
  chkxer('ZTPMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 9;
  ztpmlqt('R', 'N', 2, 2, 2, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
  chkxer('ZTPMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 11;
  ztpmlqt('R', 'N', 1, 1, 1, 1, 1, A, 1, T, 0, B, 1, C, 1, W, INFO);
  chkxer('ZTPMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 13;
  ztpmlqt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 0, C, 1, W, INFO);
  chkxer('ZTPMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 15;
  ztpmlqt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 0, W, INFO);
  chkxer('ZTPMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // Print a summary line.

  alaesm(PATH, infoc.OK.value, NOUT);
}
