// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/ztpmqrt.dart';
import 'package:lapack/src/ztpqrt.dart';
import 'package:lapack/src/ztpqrt2.dart';

import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void zerrqrtp(final String PATH, final Nout NUNIT) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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

  // Error exits for TPQRT factorization

  // ZTPQRT

  srnamc.SRNAMT = 'ZTPQRT';
  infoc.INFOT = 1;
  ztpqrt(-1, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO);
  chkxer('ZTPQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  ztpqrt(1, -1, 0, 1, A, 1, B, 1, T, 1, W, INFO);
  chkxer('ZTPQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  ztpqrt(0, 1, -1, 1, A, 1, B, 1, T, 1, W, INFO);
  chkxer('ZTPQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  ztpqrt(0, 1, 1, 1, A, 1, B, 1, T, 1, W, INFO);
  chkxer('ZTPQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  ztpqrt(0, 1, 0, 0, A, 1, B, 1, T, 1, W, INFO);
  chkxer('ZTPQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  ztpqrt(0, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO);
  chkxer('ZTPQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 6;
  ztpqrt(1, 2, 0, 2, A, 1, B, 1, T, 1, W, INFO);
  chkxer('ZTPQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  ztpqrt(2, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO);
  chkxer('ZTPQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 10;
  ztpqrt(2, 2, 1, 2, A, 2, B, 2, T, 1, W, INFO);
  chkxer('ZTPQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // ZTPQRT2

  srnamc.SRNAMT = 'ZTPQRT2';
  infoc.INFOT = 1;
  ztpqrt2(-1, 0, 0, A, 1, B, 1, T, 1, INFO);
  chkxer('ZTPQRT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  ztpqrt2(0, -1, 0, A, 1, B, 1, T, 1, INFO);
  chkxer('ZTPQRT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  ztpqrt2(0, 0, -1, A, 1, B, 1, T, 1, INFO);
  chkxer('ZTPQRT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  ztpqrt2(2, 2, 0, A, 1, B, 2, T, 2, INFO);
  chkxer('ZTPQRT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  ztpqrt2(2, 2, 0, A, 2, B, 1, T, 2, INFO);
  chkxer('ZTPQRT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 9;
  ztpqrt2(2, 2, 0, A, 2, B, 2, T, 1, INFO);
  chkxer('ZTPQRT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // ZTPMQRT

  srnamc.SRNAMT = 'ZTPMQRT';
  infoc.INFOT = 1;
  ztpmqrt('/', 'N', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
  chkxer('ZTPMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  ztpmqrt('L', '/', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
  chkxer('ZTPMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  ztpmqrt('L', 'N', -1, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
  chkxer('ZTPMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  ztpmqrt('L', 'N', 0, -1, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
  chkxer('ZTPMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  ztpmqrt('L', 'N', 0, 0, -1, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
  infoc.INFOT = 6;
  ztpmqrt('L', 'N', 0, 0, 0, -1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
  chkxer('ZTPMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  ztpmqrt('L', 'N', 0, 0, 0, 0, 0, A, 1, T, 1, B, 1, C, 1, W, INFO);
  chkxer('ZTPMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 9;
  ztpmqrt('R', 'N', 1, 2, 1, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
  chkxer('ZTPMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 9;
  ztpmqrt('L', 'N', 2, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
  chkxer('ZTPMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 11;
  ztpmqrt('R', 'N', 1, 1, 1, 1, 1, A, 1, T, 0, B, 1, C, 1, W, INFO);
  chkxer('ZTPMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 13;
  ztpmqrt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 0, C, 1, W, INFO);
  chkxer('ZTPMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 15;
  ztpmqrt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 0, W, INFO);
  chkxer('ZTPMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // Print a summary line.

  alaesm(PATH, infoc.OK.value, NOUT);
}
