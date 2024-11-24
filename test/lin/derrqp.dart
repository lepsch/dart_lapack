// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgeqp3.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void derrqp(final String PATH, final Nout NUNIT) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 3;
  final IP = Array<int>(NMAX);
  final A = Matrix<double>(NMAX, NMAX),
      TAU = Array<double>(NMAX),
      W = Array<double>(3 * NMAX + 1);
  final INFO = Box(0);

  infoc.NOUT = NUNIT;
  infoc.NOUT.println();
  final C2 = PATH.substring(1, 3);
  final LW = 3 * NMAX + 1;
  A[1][1] = 1.0;
  A[1][2] = 2.0;
  A[2][2] = 3.0;
  A[2][1] = 4.0;
  infoc.OK.value = true;

  if (lsamen(2, C2, 'QP')) {
    // Test error exits for QR factorization with pivoting

    // DGEQP3

    srnamc.SRNAMT = 'DGEQP3';
    infoc.INFOT = 1;
    dgeqp3(-1, 0, A, 1, IP, TAU, W, LW, INFO);
    chkxer('DGEQP3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgeqp3(1, -1, A, 1, IP, TAU, W, LW, INFO);
    chkxer('DGEQP3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    dgeqp3(2, 3, A, 1, IP, TAU, W, LW, INFO);
    chkxer('DGEQP3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    dgeqp3(2, 2, A, 2, IP, TAU, W, LW - 10, INFO);
    chkxer('DGEQP3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  }

  // Print a summary line.

  alaesm(PATH, infoc.OK.value, infoc.NOUT);
}
