// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/lsamen.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/zgels.dart';
import 'package:dart_lapack/src/zgelsd.dart';
import 'package:dart_lapack/src/zgelss.dart';
import 'package:dart_lapack/src/zgelst.dart';
import 'package:dart_lapack/src/zgelsy.dart';
import 'package:dart_lapack/src/zgetsls.dart';

import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void zerrls(final String PATH, final Nout NUNIT) {
  const NMAX = 2;
  final IP = Array<int>(NMAX);
  final RW = Array<double>(NMAX), S = Array<double>(NMAX);
  final A = Matrix<Complex>(NMAX, NMAX),
      B = Matrix<Complex>(NMAX, NMAX),
      W = Array<Complex>(NMAX);
  final INFO = Box(0);

  final NOUT = NUNIT;
  final C2 = PATH.substring(1, 3);
  A[1][1] = Complex(1.0, 0.0);
  A[1][2] = Complex(2.0, 0.0);
  A[2][2] = Complex(3.0, 0.0);
  A[2][1] = Complex(4.0, 0.0);
  infoc.OK.value = true;
  NOUT.println();

  // Test error exits for the least squares driver routines.

  if (lsamen(2, C2, 'LS')) {
    // ZGELS

    srnamc.SRNAMT = 'ZGELS';
    infoc.INFOT = 1;
    zgels('/', 0, 0, 0, A, 1, B, 1, W, 1, INFO);
    chkxer('ZGELS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgels('N', -1, 0, 0, A, 1, B, 1, W, 1, INFO);
    chkxer('ZGELS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgels('N', 0, -1, 0, A, 1, B, 1, W, 1, INFO);
    chkxer('ZGELS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgels('N', 0, 0, -1, A, 1, B, 1, W, 1, INFO);
    chkxer('ZGELS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zgels('N', 2, 0, 0, A, 1, B, 2, W, 2, INFO);
    chkxer('ZGELS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zgels('N', 2, 0, 0, A, 2, B, 1, W, 2, INFO);
    chkxer('ZGELS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zgels('N', 0, 2, 0, A, 1, B, 1, W, 2, INFO);
    chkxer('ZGELS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zgels('N', 1, 1, 0, A, 1, B, 1, W, 1, INFO);
    chkxer('ZGELS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGELST

    srnamc.SRNAMT = 'ZGELST';
    infoc.INFOT = 1;
    zgelst('/', 0, 0, 0, A, 1, B, 1, W, 1, INFO);
    chkxer('ZGELST', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgelst('N', -1, 0, 0, A, 1, B, 1, W, 1, INFO);
    chkxer('ZGELST', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgelst('N', 0, -1, 0, A, 1, B, 1, W, 1, INFO);
    chkxer('ZGELST', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgelst('N', 0, 0, -1, A, 1, B, 1, W, 1, INFO);
    chkxer('ZGELST', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zgelst('N', 2, 0, 0, A, 1, B, 2, W, 2, INFO);
    chkxer('ZGELST', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zgelst('N', 2, 0, 0, A, 2, B, 1, W, 2, INFO);
    chkxer('ZGELST', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zgelst('N', 0, 2, 0, A, 1, B, 1, W, 2, INFO);
    chkxer('ZGELST', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zgelst('N', 1, 1, 0, A, 1, B, 1, W, 1, INFO);
    chkxer('ZGELST', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGETSLS

    srnamc.SRNAMT = 'ZGETSLS';
    infoc.INFOT = 1;
    zgetsls('/', 0, 0, 0, A, 1, B, 1, W, 1, INFO);
    chkxer('ZGETSLS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgetsls('N', -1, 0, 0, A, 1, B, 1, W, 1, INFO);
    chkxer('ZGETSLS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgetsls('N', 0, -1, 0, A, 1, B, 1, W, 1, INFO);
    chkxer('ZGETSLS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgetsls('N', 0, 0, -1, A, 1, B, 1, W, 1, INFO);
    chkxer('ZGETSLS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zgetsls('N', 2, 0, 0, A, 1, B, 2, W, 2, INFO);
    chkxer('ZGETSLS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zgetsls('N', 2, 0, 0, A, 2, B, 1, W, 2, INFO);
    chkxer('ZGETSLS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zgetsls('N', 0, 2, 0, A, 1, B, 1, W, 2, INFO);
    chkxer('ZGETSLS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGELSS

    srnamc.SRNAMT = 'ZGELSS';
    infoc.INFOT = 1;
    final IRNK = Box(0);
    final RCOND = 0.0;
    zgelss(-1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, RW, INFO);
    chkxer('ZGELSS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgelss(0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, RW, INFO);
    chkxer('ZGELSS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgelss(0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 1, RW, INFO);
    chkxer('ZGELSS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgelss(2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 2, RW, INFO);
    chkxer('ZGELSS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zgelss(2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 2, RW, INFO);
    chkxer('ZGELSS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGELSY

    srnamc.SRNAMT = 'ZGELSY';
    infoc.INFOT = 1;
    zgelsy(-1, 0, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO);
    chkxer('ZGELSY', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgelsy(0, -1, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO);
    chkxer('ZGELSY', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgelsy(0, 0, -1, A, 1, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO);
    chkxer('ZGELSY', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgelsy(2, 0, 0, A, 1, B, 2, IP, RCOND, IRNK, W, 10, RW, INFO);
    chkxer('ZGELSY', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zgelsy(2, 0, 0, A, 2, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO);
    chkxer('ZGELSY', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    zgelsy(0, 3, 0, A, 1, B, 3, IP, RCOND, IRNK, W, 1, RW, INFO);
    chkxer('ZGELSY', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGELSD

    srnamc.SRNAMT = 'ZGELSD';
    infoc.INFOT = 1;
    zgelsd(-1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO);
    chkxer('ZGELSD', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgelsd(0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO);
    chkxer('ZGELSD', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgelsd(0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO);
    chkxer('ZGELSD', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgelsd(2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 10, RW, IP, INFO);
    chkxer('ZGELSD', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zgelsd(2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO);
    chkxer('ZGELSD', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    zgelsd(2, 2, 1, A, 2, B, 2, S, RCOND, IRNK, W, 1, RW, IP, INFO);
    chkxer('ZGELSD', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  }

  // Print a summary line.

  alaesm(PATH, infoc.OK.value, NOUT);
}
