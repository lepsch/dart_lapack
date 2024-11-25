// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/zunhr_col.dart';

import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void zerrunhr_col(String PATH, Nout NUNIT) {
  const NMAX = 2;
  final A = Matrix<Complex>(NMAX, NMAX),
      T = Matrix<Complex>(NMAX, NMAX),
      D = Array<Complex>(NMAX);
  final INFO = Box(0);

  infoc.NOUT = NUNIT;
  infoc.NOUT.println();

  // Set the variables to innocuous values.
  for (var J = 1; J <= NMAX; J++) {
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = Complex(1.0 / (I + J));
      T[I][J] = Complex(1.0 / (I + J));
    }
    D[J] = Complex(0.0, 0.0);
  }
  infoc.OK.value = true;

  // Error exits for Householder reconstruction

  // ZUNHR_COL
  srnamc.SRNAMT = 'ZUNHR_COL';

  infoc.INFOT = 1;
  zunhr_col(-1, 0, 1, A, 1, T, 1, D, INFO);
  chkxer('ZUNHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  infoc.INFOT = 2;
  zunhr_col(0, -1, 1, A, 1, T, 1, D, INFO);
  chkxer('ZUNHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  zunhr_col(1, 2, 1, A, 1, T, 1, D, INFO);
  chkxer('ZUNHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  infoc.INFOT = 3;
  zunhr_col(0, 0, -1, A, 1, T, 1, D, INFO);
  chkxer('ZUNHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  zunhr_col(0, 0, 0, A, 1, T, 1, D, INFO);
  chkxer('ZUNHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  infoc.INFOT = 5;
  zunhr_col(0, 0, 1, A, -1, T, 1, D, INFO);
  chkxer('ZUNHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  zunhr_col(0, 0, 1, A, 0, T, 1, D, INFO);
  chkxer('ZUNHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  zunhr_col(2, 0, 1, A, 1, T, 1, D, INFO);
  chkxer('ZUNHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  infoc.INFOT = 7;
  zunhr_col(0, 0, 1, A, 1, T, -1, D, INFO);
  chkxer('ZUNHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  zunhr_col(0, 0, 1, A, 1, T, 0, D, INFO);
  chkxer('ZUNHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  zunhr_col(4, 3, 2, A, 4, T, 1, D, INFO);
  chkxer('ZUNHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  // Print a summary line.
  alaesm(PATH, infoc.OK.value, infoc.NOUT);
}
