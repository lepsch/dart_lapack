// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/zgemqrt.dart';
import 'package:dart_lapack/src/zgeqrt.dart';
import 'package:dart_lapack/src/zgeqrt2.dart';
import 'package:dart_lapack/src/zgeqrt3.dart';

import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void zerrqrt(final String PATH, final Nout NUNIT) {
  const NMAX = 2;
  final A = Matrix<Complex>(NMAX, NMAX),
      T = Matrix<Complex>(NMAX, NMAX),
      W = Array<Complex>(NMAX),
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

  // Error exits for QRT factorization

  // ZGEQRT

  srnamc.SRNAMT = 'ZGEQRT';
  infoc.INFOT = 1;
  zgeqrt(-1, 0, 1, A, 1, T, 1, W, INFO);
  chkxer('ZGEQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zgeqrt(0, -1, 1, A, 1, T, 1, W, INFO);
  chkxer('ZGEQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  zgeqrt(0, 0, 0, A, 1, T, 1, W, INFO);
  chkxer('ZGEQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zgeqrt(2, 1, 1, A, 1, T, 1, W, INFO);
  chkxer('ZGEQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  zgeqrt(2, 2, 2, A, 2, T, 1, W, INFO);
  chkxer('ZGEQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // ZGEQRT2

  srnamc.SRNAMT = 'ZGEQRT2';
  infoc.INFOT = 1;
  zgeqrt2(-1, 0, A, 1, T, 1, INFO);
  chkxer('ZGEQRT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zgeqrt2(0, -1, A, 1, T, 1, INFO);
  chkxer('ZGEQRT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  zgeqrt2(2, 1, A, 1, T, 1, INFO);
  chkxer('ZGEQRT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 6;
  zgeqrt2(2, 2, A, 2, T, 1, INFO);
  chkxer('ZGEQRT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // ZGEQRT3

  srnamc.SRNAMT = 'ZGEQRT3';
  infoc.INFOT = 1;
  zgeqrt3(-1, 0, A, 1, T, 1, INFO);
  chkxer('ZGEQRT3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zgeqrt3(0, -1, A, 1, T, 1, INFO);
  chkxer('ZGEQRT3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  zgeqrt3(2, 1, A, 1, T, 1, INFO);
  chkxer('ZGEQRT3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 6;
  zgeqrt3(2, 2, A, 2, T, 1, INFO);
  chkxer('ZGEQRT3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // ZGEMQRT

  srnamc.SRNAMT = 'ZGEMQRT';
  infoc.INFOT = 1;
  zgemqrt('/', 'N', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('ZGEMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zgemqrt('L', '/', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('ZGEMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  zgemqrt('L', 'N', -1, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('ZGEMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  zgemqrt('L', 'N', 0, -1, 0, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('ZGEMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zgemqrt('L', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('ZGEMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zgemqrt('R', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('ZGEMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 6;
  zgemqrt('L', 'N', 0, 0, 0, 0, A, 1, T, 1, C, 1, W, INFO);
  chkxer('ZGEMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  zgemqrt('R', 'N', 1, 2, 1, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('ZGEMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  zgemqrt('L', 'N', 2, 1, 1, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('ZGEMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 10;
  zgemqrt('R', 'N', 1, 1, 1, 1, A, 1, T, 0, C, 1, W, INFO);
  chkxer('ZGEMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 12;
  zgemqrt('L', 'N', 1, 1, 1, 1, A, 1, T, 1, C, 0, W, INFO);
  chkxer('ZGEMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // Print a summary line.

  alaesm(PATH, infoc.OK.value, NOUT);
}
