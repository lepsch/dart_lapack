// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zgtcon.dart';
import 'package:lapack/src/zgtrfs.dart';
import 'package:lapack/src/zgttrf.dart';
import 'package:lapack/src/zgttrs.dart';
import 'package:lapack/src/zptcon.dart';
import 'package:lapack/src/zptrfs.dart';
import 'package:lapack/src/zpttrf.dart';
import 'package:lapack/src/zpttrs.dart';

import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void zerrgt(final String PATH, final Nout NUNIT) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 2;
  final IP = Array<int>(NMAX);
  final D = Array<double>(NMAX),
      DF = Array<double>(NMAX),
      R1 = Array<double>(NMAX),
      R2 = Array<double>(NMAX),
      RW = Array<double>(NMAX);
  final B = Array<Complex>(NMAX),
      DL = Array<Complex>(NMAX),
      DLF = Array<Complex>(NMAX),
      DU = Array<Complex>(NMAX),
      DU2 = Array<Complex>(NMAX),
      DUF = Array<Complex>(NMAX),
      E = Array<Complex>(NMAX),
      EF = Array<Complex>(NMAX),
      W = Array<Complex>(NMAX),
      X = Array<Complex>(NMAX);
  final INFO = Box(0);

  final NOUT = NUNIT;
  NOUT.println();
  final C2 = PATH.substring(1, 3);
  for (var I = 1; I <= NMAX; I++) {
    D[I] = 1.0;
    E[I] = 2.0.toComplex();
    DL[I] = 3.0.toComplex();
    DU[I] = 4.0.toComplex();
  }
  final ANORM = 1.0;
  infoc.OK.value = true;

  if (lsamen(2, C2, 'GT')) {
    // Test error exits for the general tridiagonal routines.

    // ZGTTRF

    srnamc.SRNAMT = 'ZGTTRF';
    infoc.INFOT = 1;
    zgttrf(-1, DL, E, DU, DU2, IP, INFO);
    chkxer('ZGTTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGTTRS

    srnamc.SRNAMT = 'ZGTTRS';
    infoc.INFOT = 1;
    zgttrs('/', 0, 0, DL, E, DU, DU2, IP, X.asMatrix(), 1, INFO);
    chkxer('ZGTTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgttrs('N', -1, 0, DL, E, DU, DU2, IP, X.asMatrix(), 1, INFO);
    chkxer('ZGTTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgttrs('N', 0, -1, DL, E, DU, DU2, IP, X.asMatrix(), 1, INFO);
    chkxer('ZGTTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zgttrs('N', 2, 1, DL, E, DU, DU2, IP, X.asMatrix(), 1, INFO);
    chkxer('ZGTTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGTRFS

    srnamc.SRNAMT = 'ZGTRFS';
    infoc.INFOT = 1;
    zgtrfs('/', 0, 0, DL, E, DU, DLF, EF, DUF, DU2, IP, B.asMatrix(), 1,
        X.asMatrix(), 1, R1, R2, W, RW, INFO);
    chkxer('ZGTRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgtrfs('N', -1, 0, DL, E, DU, DLF, EF, DUF, DU2, IP, B.asMatrix(), 1,
        X.asMatrix(), 1, R1, R2, W, RW, INFO);
    chkxer('ZGTRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgtrfs('N', 0, -1, DL, E, DU, DLF, EF, DUF, DU2, IP, B.asMatrix(), 1,
        X.asMatrix(), 1, R1, R2, W, RW, INFO);
    chkxer('ZGTRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    zgtrfs('N', 2, 1, DL, E, DU, DLF, EF, DUF, DU2, IP, B.asMatrix(), 1,
        X.asMatrix(), 2, R1, R2, W, RW, INFO);
    chkxer('ZGTRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 15;
    zgtrfs('N', 2, 1, DL, E, DU, DLF, EF, DUF, DU2, IP, B.asMatrix(), 2,
        X.asMatrix(), 1, R1, R2, W, RW, INFO);
    chkxer('ZGTRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGTCON

    srnamc.SRNAMT = 'ZGTCON';
    infoc.INFOT = 1;
    final RCOND = Box(0.0);
    zgtcon('/', 0, DL, E, DU, DU2, IP, ANORM, RCOND, W, INFO);
    chkxer('ZGTCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgtcon('I', -1, DL, E, DU, DU2, IP, ANORM, RCOND, W, INFO);
    chkxer('ZGTCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zgtcon('I', 0, DL, E, DU, DU2, IP, -ANORM, RCOND, W, INFO);
    chkxer('ZGTCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  } else if (lsamen(2, C2, 'PT')) {
    // Test error exits for the positive definite tridiagonal
    // routines.

    // ZPTTRF

    srnamc.SRNAMT = 'ZPTTRF';
    infoc.INFOT = 1;
    zpttrf(-1, D, E, INFO);
    chkxer('ZPTTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZPTTRS

    srnamc.SRNAMT = 'ZPTTRS';
    infoc.INFOT = 1;
    zpttrs('/', 1, 0, D, E, X.asMatrix(), 1, INFO);
    chkxer('ZPTTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zpttrs('U', -1, 0, D, E, X.asMatrix(), 1, INFO);
    chkxer('ZPTTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zpttrs('U', 0, -1, D, E, X.asMatrix(), 1, INFO);
    chkxer('ZPTTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zpttrs('U', 2, 1, D, E, X.asMatrix(), 1, INFO);
    chkxer('ZPTTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZPTRFS

    srnamc.SRNAMT = 'ZPTRFS';
    infoc.INFOT = 1;
    zptrfs('/', 1, 0, D, E, DF, EF, B.asMatrix(), 1, X.asMatrix(), 1, R1, R2, W,
        RW, INFO);
    chkxer('ZPTRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zptrfs('U', -1, 0, D, E, DF, EF, B.asMatrix(), 1, X.asMatrix(), 1, R1, R2,
        W, RW, INFO);
    chkxer('ZPTRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zptrfs('U', 0, -1, D, E, DF, EF, B.asMatrix(), 1, X.asMatrix(), 1, R1, R2,
        W, RW, INFO);
    chkxer('ZPTRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    zptrfs('U', 2, 1, D, E, DF, EF, B.asMatrix(), 1, X.asMatrix(), 2, R1, R2, W,
        RW, INFO);
    chkxer('ZPTRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zptrfs('U', 2, 1, D, E, DF, EF, B.asMatrix(), 2, X.asMatrix(), 1, R1, R2, W,
        RW, INFO);
    chkxer('ZPTRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZPTCON

    srnamc.SRNAMT = 'ZPTCON';
    infoc.INFOT = 1;
    final RCOND = Box(0.0);
    zptcon(-1, D, E, ANORM, RCOND, RW, INFO);
    chkxer('ZPTCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zptcon(0, D, E, -ANORM, RCOND, RW, INFO);
    chkxer('ZPTCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  }

  // Print a summary line.

  alaesm(PATH, infoc.OK.value, NOUT);
}
