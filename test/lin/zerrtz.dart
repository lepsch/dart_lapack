// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/lsamen.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/ztzrzf.dart';

import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void zerrtz(final String PATH, final Nout NUNIT) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 2;
  final A = Matrix<Complex>(NMAX, NMAX),
      TAU = Array<Complex>(NMAX),
      W = Array<Complex>(NMAX);
  final INFO = Box(0);

  final NOUT = NUNIT;
  final C2 = PATH.substring(1, 3);
  A[1][1] = Complex(1.0, -1.0);
  A[1][2] = Complex(2.0, -2.0);
  A[2][2] = Complex(3.0, -3.0);
  A[2][1] = Complex(4.0, -4.0);
  W[1] = Complex.zero;
  W[2] = Complex.zero;
  infoc.OK.value = true;

  // Test error exits for the trapezoidal routines.
  NOUT.println();
  if (lsamen(2, C2, 'TZ')) {
    // ZTZRZF

    srnamc.SRNAMT = 'ZTZRZF';
    infoc.INFOT = 1;
    ztzrzf(-1, 0, A, 1, TAU, W, 1, INFO);
    chkxer('ZTZRZF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    ztzrzf(1, 0, A, 1, TAU, W, 1, INFO);
    chkxer('ZTZRZF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    ztzrzf(2, 2, A, 1, TAU, W, 1, INFO);
    chkxer('ZTZRZF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    ztzrzf(2, 2, A, 2, TAU, W, 0, INFO);
    chkxer('ZTZRZF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    ztzrzf(2, 3, A, 2, TAU, W, 1, INFO);
    chkxer('ZTZRZF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  }

  // Print a summary line.

  alaesm(PATH, infoc.OK.value, NOUT);
}
