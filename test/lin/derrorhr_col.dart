// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';

import '../test_driver.dart';
import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void derrorhr_col(final String PATH, final Nout NUNIT, final TestDriver test) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 2;
  final A = Matrix<double>(NMAX, NMAX),
      T = Matrix<double>(NMAX, NMAX),
      D = Array<double>(NMAX);
  final INFO = Box(0);

  final NOUT = infoc.NOUT = NUNIT;
  final OK = infoc.OK;
  final LERR = infoc.LERR;
  NOUT.println();

  // Set the variables to innocuous values.

  for (var J = 1; J <= NMAX; J++) {
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = 1.0 / (I + J);
      T[I][J] = 1.0 / (I + J);
    }
    D[J] = 0.0;
  }
  OK.value = true;

  // Error exits for Householder reconstruction

  test('DORHR_COL', () {
    srnamc.SRNAMT = 'DORHR_COL';

    infoc.INFOT = 1;
    dorhr_col(-1, 0, 1, A, 1, T, 1, D, INFO);
    chkxer('DORHR_COL', infoc.INFOT, NOUT, LERR, OK, test);

    infoc.INFOT = 2;
    dorhr_col(0, -1, 1, A, 1, T, 1, D, INFO);
    chkxer('DORHR_COL', infoc.INFOT, NOUT, LERR, OK, test);
    dorhr_col(1, 2, 1, A, 1, T, 1, D, INFO);
    chkxer('DORHR_COL', infoc.INFOT, NOUT, LERR, OK, test);

    infoc.INFOT = 3;
    dorhr_col(0, 0, -1, A, 1, T, 1, D, INFO);
    chkxer('DORHR_COL', infoc.INFOT, NOUT, LERR, OK, test);

    dorhr_col(0, 0, 0, A, 1, T, 1, D, INFO);
    chkxer('DORHR_COL', infoc.INFOT, NOUT, LERR, OK, test);

    infoc.INFOT = 5;
    dorhr_col(0, 0, 1, A, -1, T, 1, D, INFO);
    chkxer('DORHR_COL', infoc.INFOT, NOUT, LERR, OK, test);

    dorhr_col(0, 0, 1, A, 0, T, 1, D, INFO);
    chkxer('DORHR_COL', infoc.INFOT, NOUT, LERR, OK, test);

    dorhr_col(2, 0, 1, A, 1, T, 1, D, INFO);
    chkxer('DORHR_COL', infoc.INFOT, NOUT, LERR, OK, test);

    infoc.INFOT = 7;
    dorhr_col(0, 0, 1, A, 1, T, -1, D, INFO);
    chkxer('DORHR_COL', infoc.INFOT, NOUT, LERR, OK, test);

    dorhr_col(0, 0, 1, A, 1, T, 0, D, INFO);
    chkxer('DORHR_COL', infoc.INFOT, NOUT, LERR, OK, test);

    dorhr_col(4, 3, 2, A, 4, T, 1, D, INFO);
    chkxer('DORHR_COL', infoc.INFOT, NOUT, LERR, OK, test);
  });

  // Print a summary line.
  alaesm(PATH, OK.value, NOUT);
}
