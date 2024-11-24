// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dsgesv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import 'chkxer.dart';
import 'common.dart';

void derrab(final Nout NUNIT) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 4;
  final IP = Array<int>(NMAX);
  final A = Matrix<double>(NMAX, NMAX),
      AF = Matrix<double>(NMAX, NMAX),
      B = Array<double>(NMAX),
      C = Array<double>(NMAX),
      R = Array<double>(NMAX),
      R1 = Array<double>(NMAX),
      R2 = Array<double>(NMAX),
      W = Array<double>(2 * NMAX),
      X = Array<double>(NMAX);
  final WORK = Array<double>(1);
  final SWORK = Array<double>(1);
  final INFO = Box(0), ITER = Box(0);

  infoc.NOUT = NUNIT;
  infoc.NOUT.println();

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
    C[J] = 0.0;
    R[J] = 0.0;
    IP[J] = J;
  }
  infoc.OK.value = true;

  srnamc.SRNAMT = 'DSGESV';
  infoc.INFOT = 1;
  dsgesv(-1, 0, A, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1, WORK.asMatrix(),
      SWORK, ITER, INFO);
  chkxer('DSGESV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dsgesv(0, -1, A, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1, WORK.asMatrix(),
      SWORK, ITER, INFO);
  chkxer('DSGESV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  dsgesv(2, 1, A, 1, IP, B.asMatrix(), 2, X.asMatrix(), 2, WORK.asMatrix(),
      SWORK, ITER, INFO);
  chkxer('DSGESV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  dsgesv(2, 1, A, 2, IP, B.asMatrix(), 1, X.asMatrix(), 2, WORK.asMatrix(),
      SWORK, ITER, INFO);
  chkxer('DSGESV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 9;
  dsgesv(2, 1, A, 2, IP, B.asMatrix(), 2, X.asMatrix(), 1, WORK.asMatrix(),
      SWORK, ITER, INFO);
  chkxer('DSGESV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  // Print a summary line.

  if (infoc.OK.value) {
    infoc.NOUT.println(' DSGESV drivers passed the tests of the error exits');
  } else {
    infoc.NOUT
        .println(' *** DSGESV drivers failed the tests of the error exits ***');
  }
}
