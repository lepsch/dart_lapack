// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dtzrzf.dart';
import 'package:dart_lapack/src/lsamen.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';

import '../test_driver.dart';
import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void derrtz(final String PATH, final Nout NUNIT, final TestDriver test) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 2;
  final A = Matrix<double>(NMAX, NMAX),
      TAU = Array<double>(NMAX),
      W = Array<double>(NMAX);
  final INFO = Box(0);

  infoc.NOUT = NUNIT;
  infoc.NOUT.println();
  final C2 = PATH.substring(1, 3);
  A[1][1] = 1.0;
  A[1][2] = 2.0;
  A[2][2] = 3.0;
  A[2][1] = 4.0;
  W[1] = 0.0;
  W[2] = 0.0;
  infoc.OK.value = true;

  if (lsamen(2, C2, 'TZ')) {
    // Test error exits for the trapezoidal routines.

    test('DTZRZF', () {
      srnamc.SRNAMT = 'DTZRZF';
      infoc.INFOT = 1;
      dtzrzf(-1, 0, A, 1, TAU, W, 1, INFO);
      chkxer('DTZRZF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dtzrzf(1, 0, A, 1, TAU, W, 1, INFO);
      chkxer('DTZRZF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dtzrzf(2, 2, A, 1, TAU, W, 1, INFO);
      chkxer('DTZRZF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dtzrzf(2, 2, A, 2, TAU, W, 0, INFO);
      chkxer('DTZRZF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dtzrzf(2, 3, A, 2, TAU, W, 1, INFO);
      chkxer('DTZRZF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    });
  }

  // Print a summary line.
  alaesm(PATH, infoc.OK.value, infoc.NOUT);
}
