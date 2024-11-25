// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/lsamen.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/zgeqp3.dart';

import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void zerrqp(final String PATH, final Nout NUNIT) {
  const NMAX = 3;
  final IP = Array<int>(NMAX);
  final RW = Array<double>(2 * NMAX);
  final A = Matrix<Complex>(NMAX, NMAX),
      TAU = Array<Complex>(NMAX),
      W = Array<Complex>(2 * NMAX + 3 * NMAX);
  final INFO = Box(0);

  final NOUT = NUNIT;
  final C2 = PATH.substring(1, 3);
  final LW = NMAX + 1;
  A[1][1] = Complex(1.0, -1.0);
  A[1][2] = Complex(2.0, -2.0);
  A[2][2] = Complex(3.0, -3.0);
  A[2][1] = Complex(4.0, -4.0);
  infoc.OK.value = true;
  NOUT.println();

  // Test error exits for QR factorization with pivoting

  if (lsamen(2, C2, 'QP')) {
    // ZGEQP3

    srnamc.SRNAMT = 'ZGEQP3';
    infoc.INFOT = 1;
    zgeqp3(-1, 0, A, 1, IP, TAU, W, LW, RW, INFO);
    chkxer('ZGEQP3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgeqp3(1, -1, A, 1, IP, TAU, W, LW, RW, INFO);
    chkxer('ZGEQP3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgeqp3(2, 3, A, 1, IP, TAU, W, LW, RW, INFO);
    chkxer('ZGEQP3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zgeqp3(2, 2, A, 2, IP, TAU, W, LW - 10, RW, INFO);
    chkxer('ZGEQP3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  }

  // Print a summary line.

  alaesm(PATH, infoc.OK.value, NOUT);
}
