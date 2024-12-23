// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/lsamen.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/zgbsv.dart';
import 'package:dart_lapack/src/zgbsvx.dart';
import 'package:dart_lapack/src/zgesv.dart';
import 'package:dart_lapack/src/zgesvx.dart';
import 'package:dart_lapack/src/zgtsv.dart';
import 'package:dart_lapack/src/zgtsvx.dart';
import 'package:dart_lapack/src/zhesv.dart';
import 'package:dart_lapack/src/zhesv_aa.dart';
import 'package:dart_lapack/src/zhesv_aa_2stage.dart';
import 'package:dart_lapack/src/zhesv_rk.dart';
import 'package:dart_lapack/src/zhesv_rook.dart';
import 'package:dart_lapack/src/zhesvx.dart';
import 'package:dart_lapack/src/zhpsv.dart';
import 'package:dart_lapack/src/zhpsvx.dart';
import 'package:dart_lapack/src/zpbsv.dart';
import 'package:dart_lapack/src/zpbsvx.dart';
import 'package:dart_lapack/src/zposv.dart';
import 'package:dart_lapack/src/zposvx.dart';
import 'package:dart_lapack/src/zppsv.dart';
import 'package:dart_lapack/src/zppsvx.dart';
import 'package:dart_lapack/src/zptsv.dart';
import 'package:dart_lapack/src/zptsvx.dart';
import 'package:dart_lapack/src/zspsv.dart';
import 'package:dart_lapack/src/zspsvx.dart';
import 'package:dart_lapack/src/zsysv.dart';
import 'package:dart_lapack/src/zsysv_aa.dart';
import 'package:dart_lapack/src/zsysv_aa_2stage.dart';
import 'package:dart_lapack/src/zsysv_rk.dart';
import 'package:dart_lapack/src/zsysv_rook.dart';
import 'package:dart_lapack/src/zsysvx.dart';

import 'chkxer.dart';
import 'common.dart';

void zerrvx(final String PATH, final Nout NUNIT) {
  const NMAX = 4;
  String EQ;
  final IP = Array<int>(NMAX);
  final C = Array<double>(NMAX),
      R = Array<double>(NMAX),
      R1 = Array<double>(NMAX),
      R2 = Array<double>(NMAX),
      RF = Array<double>(NMAX),
      RW = Array<double>(NMAX);
  final A = Matrix<Complex>(NMAX, NMAX),
      AF = Matrix<Complex>(NMAX, NMAX),
      B = Array<Complex>(NMAX),
      E = Array<Complex>(NMAX),
      W = Array<Complex>(2 * NMAX),
      X = Array<Complex>(NMAX);
  final INFO = Box(0);

  final NOUT = NUNIT;
  NOUT.println();
  final C2 = PATH.substring(1, 3);

  // Set the variables to innocuous values.

  for (var J = 1; J <= NMAX; J++) {
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = Complex(1.0 / (I + J), -1.0 / (I + J));
      AF[I][J] = Complex(1.0 / (I + J), -1.0 / (I + J));
    }
    B[J] = Complex.zero;
    E[J] = Complex.zero;
    R1[J] = 0.0;
    R2[J] = 0.0;
    W[J] = Complex.zero;
    X[J] = Complex.zero;
    C[J] = 0.0;
    R[J] = 0.0;
    IP[J] = J;
  }
  EQ = ' ';
  infoc.OK.value = true;

  if (lsamen(2, C2, 'GE')) {
    // ZGESV

    srnamc.SRNAMT = 'ZGESV';
    infoc.INFOT = 1;
    zgesv(-1, 0, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('ZGESV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgesv(0, -1, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('ZGESV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgesv(2, 1, A, 1, IP, B.asMatrix(), 2, INFO);
    chkxer('ZGESV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zgesv(2, 1, A, 2, IP, B.asMatrix(), 1, INFO);
    chkxer('ZGESV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGESVX

    srnamc.SRNAMT = 'ZGESVX';
    infoc.INFOT = 1;
    final RCOND = Box(0.0);
    zgesvx('/', 'N', 0, 0, A, 1, AF, 1, IP, Box(EQ), R, C, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGESVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgesvx('N', '/', 0, 0, A, 1, AF, 1, IP, Box(EQ), R, C, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGESVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgesvx('N', 'N', -1, 0, A, 1, AF, 1, IP, Box(EQ), R, C, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGESVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgesvx('N', 'N', 0, -1, A, 1, AF, 1, IP, Box(EQ), R, C, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGESVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zgesvx('N', 'N', 2, 1, A, 1, AF, 2, IP, Box(EQ), R, C, B.asMatrix(), 2,
        X.asMatrix(), 2, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGESVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zgesvx('N', 'N', 2, 1, A, 2, AF, 1, IP, Box(EQ), R, C, B.asMatrix(), 2,
        X.asMatrix(), 2, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGESVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    EQ = '/';
    zgesvx('F', 'N', 0, 0, A, 1, AF, 1, IP, Box(EQ), R, C, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGESVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    EQ = 'R';
    zgesvx('F', 'N', 1, 0, A, 1, AF, 1, IP, Box(EQ), R, C, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGESVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    EQ = 'C';
    zgesvx('F', 'N', 1, 0, A, 1, AF, 1, IP, Box(EQ), R, C, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGESVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 14;
    zgesvx('N', 'N', 2, 1, A, 2, AF, 2, IP, Box(EQ), R, C, B.asMatrix(), 1,
        X.asMatrix(), 2, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGESVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 16;
    zgesvx('N', 'N', 2, 1, A, 2, AF, 2, IP, Box(EQ), R, C, B.asMatrix(), 2,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGESVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  } else if (lsamen(2, C2, 'GB')) {
    // ZGBSV

    srnamc.SRNAMT = 'ZGBSV';
    infoc.INFOT = 1;
    zgbsv(-1, 0, 0, 0, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('ZGBSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgbsv(1, -1, 0, 0, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('ZGBSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgbsv(1, 0, -1, 0, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('ZGBSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgbsv(0, 0, 0, -1, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('ZGBSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zgbsv(1, 1, 1, 0, A, 3, IP, B.asMatrix(), 1, INFO);
    chkxer('ZGBSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    zgbsv(2, 0, 0, 0, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('ZGBSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGBSVX

    srnamc.SRNAMT = 'ZGBSVX';
    infoc.INFOT = 1;
    final RCOND = Box(0.0);
    zgbsvx('/', 'N', 0, 0, 0, 0, A, 1, AF, 1, IP, Box(EQ), R, C, B.asMatrix(),
        1, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgbsvx('N', '/', 0, 0, 0, 0, A, 1, AF, 1, IP, Box(EQ), R, C, B.asMatrix(),
        1, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgbsvx('N', 'N', -1, 0, 0, 0, A, 1, AF, 1, IP, Box(EQ), R, C, B.asMatrix(),
        1, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgbsvx('N', 'N', 1, -1, 0, 0, A, 1, AF, 1, IP, Box(EQ), R, C, B.asMatrix(),
        1, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgbsvx('N', 'N', 1, 0, -1, 0, A, 1, AF, 1, IP, Box(EQ), R, C, B.asMatrix(),
        1, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zgbsvx('N', 'N', 0, 0, 0, -1, A, 1, AF, 1, IP, Box(EQ), R, C, B.asMatrix(),
        1, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zgbsvx('N', 'N', 1, 1, 1, 0, A, 2, AF, 4, IP, Box(EQ), R, C, B.asMatrix(),
        1, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zgbsvx('N', 'N', 1, 1, 1, 0, A, 3, AF, 3, IP, Box(EQ), R, C, B.asMatrix(),
        1, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    EQ = '/';
    zgbsvx('F', 'N', 0, 0, 0, 0, A, 1, AF, 1, IP, Box(EQ), R, C, B.asMatrix(),
        1, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    EQ = 'R';
    zgbsvx('F', 'N', 1, 0, 0, 0, A, 1, AF, 1, IP, Box(EQ), R, C, B.asMatrix(),
        1, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 14;
    EQ = 'C';
    zgbsvx('F', 'N', 1, 0, 0, 0, A, 1, AF, 1, IP, Box(EQ), R, C, B.asMatrix(),
        1, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 16;
    zgbsvx('N', 'N', 2, 0, 0, 0, A, 1, AF, 1, IP, Box(EQ), R, C, B.asMatrix(),
        1, X.asMatrix(), 2, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 18;
    zgbsvx('N', 'N', 2, 0, 0, 0, A, 1, AF, 1, IP, Box(EQ), R, C, B.asMatrix(),
        2, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZGBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  } else if (lsamen(2, C2, 'GT')) {
    // ZGTSV

    srnamc.SRNAMT = 'ZGTSV';
    infoc.INFOT = 1;
    zgtsv(-1, 0, A(1, 1).asArray(), A(1, 2).asArray(), A(1, 3).asArray(),
        B.asMatrix(), 1, INFO);
    chkxer('ZGTSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgtsv(0, -1, A(1, 1).asArray(), A(1, 2).asArray(), A(1, 3).asArray(),
        B.asMatrix(), 1, INFO);
    chkxer('ZGTSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zgtsv(2, 0, A(1, 1).asArray(), A(1, 2).asArray(), A(1, 3).asArray(),
        B.asMatrix(), 1, INFO);
    chkxer('ZGTSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGTSVX

    srnamc.SRNAMT = 'ZGTSVX';
    infoc.INFOT = 1;
    final RCOND = Box(0.0);
    zgtsvx(
        '/',
        'N',
        0,
        0,
        A(1, 1).asArray(),
        A(1, 2).asArray(),
        A(1, 3).asArray(),
        AF(1, 1).asArray(),
        AF(1, 2).asArray(),
        AF(1, 3).asArray(),
        AF(1, 4).asArray(),
        IP,
        B.asMatrix(),
        1,
        X.asMatrix(),
        1,
        RCOND,
        R1,
        R2,
        W,
        RW,
        INFO);
    chkxer('ZGTSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgtsvx(
        'N',
        '/',
        0,
        0,
        A(1, 1).asArray(),
        A(1, 2).asArray(),
        A(1, 3).asArray(),
        AF(1, 1).asArray(),
        AF(1, 2).asArray(),
        AF(1, 3).asArray(),
        AF(1, 4).asArray(),
        IP,
        B.asMatrix(),
        1,
        X.asMatrix(),
        1,
        RCOND,
        R1,
        R2,
        W,
        RW,
        INFO);
    chkxer('ZGTSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgtsvx(
        'N',
        'N',
        -1,
        0,
        A(1, 1).asArray(),
        A(1, 2).asArray(),
        A(1, 3).asArray(),
        AF(1, 1).asArray(),
        AF(1, 2).asArray(),
        AF(1, 3).asArray(),
        AF(1, 4).asArray(),
        IP,
        B.asMatrix(),
        1,
        X.asMatrix(),
        1,
        RCOND,
        R1,
        R2,
        W,
        RW,
        INFO);
    chkxer('ZGTSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgtsvx(
        'N',
        'N',
        0,
        -1,
        A(1, 1).asArray(),
        A(1, 2).asArray(),
        A(1, 3).asArray(),
        AF(1, 1).asArray(),
        AF(1, 2).asArray(),
        AF(1, 3).asArray(),
        AF(1, 4).asArray(),
        IP,
        B.asMatrix(),
        1,
        X.asMatrix(),
        1,
        RCOND,
        R1,
        R2,
        W,
        RW,
        INFO);
    chkxer('ZGTSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 14;
    zgtsvx(
        'N',
        'N',
        2,
        0,
        A(1, 1).asArray(),
        A(1, 2).asArray(),
        A(1, 3).asArray(),
        AF(1, 1).asArray(),
        AF(1, 2).asArray(),
        AF(1, 3).asArray(),
        AF(1, 4).asArray(),
        IP,
        B.asMatrix(),
        1,
        X.asMatrix(),
        2,
        RCOND,
        R1,
        R2,
        W,
        RW,
        INFO);
    chkxer('ZGTSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 16;
    zgtsvx(
        'N',
        'N',
        2,
        0,
        A(1, 1).asArray(),
        A(1, 2).asArray(),
        A(1, 3).asArray(),
        AF(1, 1).asArray(),
        AF(1, 2).asArray(),
        AF(1, 3).asArray(),
        AF(1, 4).asArray(),
        IP,
        B.asMatrix(),
        2,
        X.asMatrix(),
        1,
        RCOND,
        R1,
        R2,
        W,
        RW,
        INFO);
    chkxer('ZGTSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  } else if (lsamen(2, C2, 'PO')) {
    // ZPOSV

    srnamc.SRNAMT = 'ZPOSV';
    infoc.INFOT = 1;
    zposv('/', 0, 0, A, 1, B.asMatrix(), 1, INFO);
    chkxer('ZPOSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zposv('U', -1, 0, A, 1, B.asMatrix(), 1, INFO);
    chkxer('ZPOSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zposv('U', 0, -1, A, 1, B.asMatrix(), 1, INFO);
    chkxer('ZPOSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zposv('U', 2, 0, A, 1, B.asMatrix(), 2, INFO);
    chkxer('ZPOSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zposv('U', 2, 0, A, 2, B.asMatrix(), 1, INFO);
    chkxer('ZPOSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZPOSVX

    srnamc.SRNAMT = 'ZPOSVX';
    infoc.INFOT = 1;
    final RCOND = Box(0.0);
    zposvx('/', 'U', 0, 0, A, 1, AF, 1, Box(EQ), C, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPOSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zposvx('N', '/', 0, 0, A, 1, AF, 1, Box(EQ), C, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPOSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zposvx('N', 'U', -1, 0, A, 1, AF, 1, Box(EQ), C, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPOSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zposvx('N', 'U', 0, -1, A, 1, AF, 1, Box(EQ), C, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPOSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zposvx('N', 'U', 2, 0, A, 1, AF, 2, Box(EQ), C, B.asMatrix(), 2,
        X.asMatrix(), 2, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPOSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zposvx('N', 'U', 2, 0, A, 2, AF, 1, Box(EQ), C, B.asMatrix(), 2,
        X.asMatrix(), 2, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPOSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    EQ = '/';
    zposvx('F', 'U', 0, 0, A, 1, AF, 1, Box(EQ), C, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPOSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    EQ = 'Y';
    zposvx('F', 'U', 1, 0, A, 1, AF, 1, Box(EQ), C, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPOSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    zposvx('N', 'U', 2, 0, A, 2, AF, 2, Box(EQ), C, B.asMatrix(), 1,
        X.asMatrix(), 2, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPOSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 14;
    zposvx('N', 'U', 2, 0, A, 2, AF, 2, Box(EQ), C, B.asMatrix(), 2,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPOSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  } else if (lsamen(2, C2, 'PP')) {
    // ZPPSV

    srnamc.SRNAMT = 'ZPPSV';
    infoc.INFOT = 1;
    zppsv('/', 0, 0, A.asArray(), B.asMatrix(), 1, INFO);
    chkxer('ZPPSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zppsv('U', -1, 0, A.asArray(), B.asMatrix(), 1, INFO);
    chkxer('ZPPSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zppsv('U', 0, -1, A.asArray(), B.asMatrix(), 1, INFO);
    chkxer('ZPPSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zppsv('U', 2, 0, A.asArray(), B.asMatrix(), 1, INFO);
    chkxer('ZPPSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZPPSVX

    srnamc.SRNAMT = 'ZPPSVX';
    infoc.INFOT = 1;
    final RCOND = Box(0.0);
    zppsvx('/', 'U', 0, 0, A.asArray(), AF.asArray(), Box(EQ), C, B.asMatrix(),
        1, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPPSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zppsvx('N', '/', 0, 0, A.asArray(), AF.asArray(), Box(EQ), C, B.asMatrix(),
        1, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPPSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zppsvx('N', 'U', -1, 0, A.asArray(), AF.asArray(), Box(EQ), C, B.asMatrix(),
        1, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPPSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zppsvx('N', 'U', 0, -1, A.asArray(), AF.asArray(), Box(EQ), C, B.asMatrix(),
        1, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPPSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    EQ = '/';
    zppsvx('F', 'U', 0, 0, A.asArray(), AF.asArray(), Box(EQ), C, B.asMatrix(),
        1, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPPSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    EQ = 'Y';
    zppsvx('F', 'U', 1, 0, A.asArray(), AF.asArray(), Box(EQ), C, B.asMatrix(),
        1, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPPSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zppsvx('N', 'U', 2, 0, A.asArray(), AF.asArray(), Box(EQ), C, B.asMatrix(),
        1, X.asMatrix(), 2, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPPSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    zppsvx('N', 'U', 2, 0, A.asArray(), AF.asArray(), Box(EQ), C, B.asMatrix(),
        2, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPPSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  } else if (lsamen(2, C2, 'PB')) {
    // ZPBSV

    srnamc.SRNAMT = 'ZPBSV';
    infoc.INFOT = 1;
    zpbsv('/', 0, 0, 0, A, 1, B.asMatrix(), 1, INFO);
    chkxer('ZPBSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zpbsv('U', -1, 0, 0, A, 1, B.asMatrix(), 1, INFO);
    chkxer('ZPBSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zpbsv('U', 1, -1, 0, A, 1, B.asMatrix(), 1, INFO);
    chkxer('ZPBSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zpbsv('U', 0, 0, -1, A, 1, B.asMatrix(), 1, INFO);
    chkxer('ZPBSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zpbsv('U', 1, 1, 0, A, 1, B.asMatrix(), 2, INFO);
    chkxer('ZPBSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zpbsv('U', 2, 0, 0, A, 1, B.asMatrix(), 1, INFO);
    chkxer('ZPBSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZPBSVX

    srnamc.SRNAMT = 'ZPBSVX';
    infoc.INFOT = 1;
    final RCOND = Box(0.0);
    zpbsvx('/', 'U', 0, 0, 0, A, 1, AF, 1, Box(EQ), C, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zpbsvx('N', '/', 0, 0, 0, A, 1, AF, 1, Box(EQ), C, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zpbsvx('N', 'U', -1, 0, 0, A, 1, AF, 1, Box(EQ), C, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zpbsvx('N', 'U', 1, -1, 0, A, 1, AF, 1, Box(EQ), C, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zpbsvx('N', 'U', 0, 0, -1, A, 1, AF, 1, Box(EQ), C, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zpbsvx('N', 'U', 1, 1, 0, A, 1, AF, 2, Box(EQ), C, B.asMatrix(), 2,
        X.asMatrix(), 2, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    zpbsvx('N', 'U', 1, 1, 0, A, 2, AF, 1, Box(EQ), C, B.asMatrix(), 2,
        X.asMatrix(), 2, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    EQ = '/';
    zpbsvx('F', 'U', 0, 0, 0, A, 1, AF, 1, Box(EQ), C, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    EQ = 'Y';
    zpbsvx('F', 'U', 1, 0, 0, A, 1, AF, 1, Box(EQ), C, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    zpbsvx('N', 'U', 2, 0, 0, A, 1, AF, 1, Box(EQ), C, B.asMatrix(), 1,
        X.asMatrix(), 2, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 15;
    zpbsvx('N', 'U', 2, 0, 0, A, 1, AF, 1, Box(EQ), C, B.asMatrix(), 2,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPBSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  } else if (lsamen(2, C2, 'PT')) {
    // ZPTSV

    srnamc.SRNAMT = 'ZPTSV';
    infoc.INFOT = 1;
    zptsv(-1, 0, R, A(1, 1).asArray(), B.asMatrix(), 1, INFO);
    chkxer('ZPTSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zptsv(0, -1, R, A(1, 1).asArray(), B.asMatrix(), 1, INFO);
    chkxer('ZPTSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zptsv(2, 0, R, A(1, 1).asArray(), B.asMatrix(), 1, INFO);
    chkxer('ZPTSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZPTSVX

    srnamc.SRNAMT = 'ZPTSVX';
    infoc.INFOT = 1;
    final RCOND = Box(0.0);
    zptsvx('/', 0, 0, R, A(1, 1).asArray(), RF, AF(1, 1).asArray(),
        B.asMatrix(), 1, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPTSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zptsvx('N', -1, 0, R, A(1, 1).asArray(), RF, AF(1, 1).asArray(),
        B.asMatrix(), 1, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPTSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zptsvx('N', 0, -1, R, A(1, 1).asArray(), RF, AF(1, 1).asArray(),
        B.asMatrix(), 1, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPTSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    zptsvx('N', 2, 0, R, A(1, 1).asArray(), RF, AF(1, 1).asArray(),
        B.asMatrix(), 1, X.asMatrix(), 2, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPTSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zptsvx('N', 2, 0, R, A(1, 1).asArray(), RF, AF(1, 1).asArray(),
        B.asMatrix(), 2, X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZPTSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  } else if (lsamen(2, C2, 'HE')) {
    // ZHESV

    srnamc.SRNAMT = 'ZHESV';
    infoc.INFOT = 1;
    zhesv('/', 0, 0, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZHESV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zhesv('U', -1, 0, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZHESV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zhesv('U', 0, -1, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZHESV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zhesv('U', 2, 0, A, 1, IP, B.asMatrix(), 2, W, 1, INFO);
    chkxer('ZHESV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zhesv('U', 2, 0, A, 2, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZHESV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zhesv('U', 0, 0, A, 1, IP, B.asMatrix(), 1, W, 0, INFO);
    chkxer('ZHESV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zhesv('U', 0, 0, A, 1, IP, B.asMatrix(), 1, W, -2, INFO);
    chkxer('ZHESV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZHESVX

    srnamc.SRNAMT = 'ZHESVX';
    infoc.INFOT = 1;
    final RCOND = Box(0.0);
    zhesvx('/', 'U', 0, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
        RCOND, R1, R2, W, 1, RW, INFO);
    chkxer('ZHESVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zhesvx('N', '/', 0, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
        RCOND, R1, R2, W, 1, RW, INFO);
    chkxer('ZHESVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zhesvx('N', 'U', -1, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
        RCOND, R1, R2, W, 1, RW, INFO);
    chkxer('ZHESVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zhesvx('N', 'U', 0, -1, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
        RCOND, R1, R2, W, 1, RW, INFO);
    chkxer('ZHESVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zhesvx('N', 'U', 2, 0, A, 1, AF, 2, IP, B.asMatrix(), 2, X.asMatrix(), 2,
        RCOND, R1, R2, W, 4, RW, INFO);
    chkxer('ZHESVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zhesvx('N', 'U', 2, 0, A, 2, AF, 1, IP, B.asMatrix(), 2, X.asMatrix(), 2,
        RCOND, R1, R2, W, 4, RW, INFO);
    chkxer('ZHESVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zhesvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B.asMatrix(), 1, X.asMatrix(), 2,
        RCOND, R1, R2, W, 4, RW, INFO);
    chkxer('ZHESVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    zhesvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B.asMatrix(), 2, X.asMatrix(), 1,
        RCOND, R1, R2, W, 4, RW, INFO);
    chkxer('ZHESVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 18;
    zhesvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B.asMatrix(), 2, X.asMatrix(), 2,
        RCOND, R1, R2, W, 3, RW, INFO);
    chkxer('ZHESVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  } else if (lsamen(2, C2, 'HR')) {
    // ZHESV_ROOK

    srnamc.SRNAMT = 'ZHESV_ROOK';
    infoc.INFOT = 1;
    zhesv_rook('/', 0, 0, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZHESV_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zhesv_rook('U', -1, 0, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZHESV_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zhesv_rook('U', 0, -1, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZHESV_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zhesv_rook('U', 2, 0, A, 1, IP, B.asMatrix(), 2, W, 1, INFO);
    chkxer('ZHESV_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zhesv_rook('U', 2, 0, A, 2, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZHESV_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zhesv_rook('U', 0, 0, A, 1, IP, B.asMatrix(), 1, W, 0, INFO);
    chkxer('ZHESV_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zhesv_rook('U', 0, 0, A, 1, IP, B.asMatrix(), 1, W, -2, INFO);
    chkxer('ZHESV_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  } else if (lsamen(2, C2, 'HK')) {
    // ZSYSV_RK

    // Test error exits of the driver that uses factorization
    // of a Hermitian indefinite matrix with rook
    // (bounded Bunch-Kaufman) pivoting with the new storage
    // format for factors L ( or U) and D.

    // L (or U) is stored in A, diagonal of D is stored on the
    // diagonal of A, subdiagonal of D is stored in a separate array E.

    srnamc.SRNAMT = 'ZHESV_RK';
    infoc.INFOT = 1;
    zhesv_rk('/', 0, 0, A, 1, E, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZHESV_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zhesv_rk('U', -1, 0, A, 1, E, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZHESV_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zhesv_rk('U', 0, -1, A, 1, E, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZHESV_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zhesv_rk('U', 2, 0, A, 1, E, IP, B.asMatrix(), 2, W, 1, INFO);
    chkxer('ZHESV_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    zhesv_rk('U', 2, 0, A, 2, E, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZHESV_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zhesv_rk('U', 0, 0, A, 1, E, IP, B.asMatrix(), 1, W, 0, INFO);
    chkxer('ZHESV_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zhesv_rk('U', 0, 0, A, 1, E, IP, B.asMatrix(), 1, W, -2, INFO);
    chkxer('ZHESV_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  } else if (lsamen(2, C2, 'HA')) {
    // ZHESV_AASEN

    srnamc.SRNAMT = 'ZHESV_AA';
    infoc.INFOT = 1;
    zhesv_aa('/', 0, 0, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZHESV_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zhesv_aa('U', -1, 0, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZHESV_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zhesv_aa('U', 0, -1, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZHESV_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zhesv_aa('U', 2, 0, A, 1, IP, B.asMatrix(), 2, W, 1, INFO);
    chkxer('ZHESV_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zhesv_aa('U', 2, 0, A, 2, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZHESV_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zhesv_aa('U', 3, 1, A, 3, IP, B.asMatrix(), 3, W, 6, INFO);
    chkxer('ZHESV_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  } else if (lsamen(2, C2, 'H2')) {
    // ZHESV_AASEN_2STAGE

    srnamc.SRNAMT = 'ZHESV_AA_2STAGE';
    infoc.INFOT = 1;
    zhesv_aa_2stage(
        '/', 0, 0, A, 1, A.asArray(), 1, IP, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZHESV_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zhesv_aa_2stage(
        'U', -1, 0, A, 1, A.asArray(), 1, IP, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZHESV_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zhesv_aa_2stage(
        'U', 0, -1, A, 1, A.asArray(), 1, IP, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZHESV_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zhesv_aa_2stage(
        'U', 2, 1, A, 1, A.asArray(), 1, IP, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZHESV_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zhesv_aa_2stage(
        'U', 2, 1, A, 2, A.asArray(), 1, IP, IP, B.asMatrix(), 2, W, 1, INFO);
    chkxer('ZHESV_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zhesv_aa_2stage(
        'U', 2, 1, A, 2, A.asArray(), 8, IP, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZHESV_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    zhesv_aa_2stage(
        'U', 2, 1, A, 2, A.asArray(), 8, IP, IP, B.asMatrix(), 2, W, 1, INFO);
    chkxer('ZHESV_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  } else if (lsamen(2, C2, 'SA')) {
    // ZSYSV_AASEN

    srnamc.SRNAMT = 'ZSYSV_AA';
    infoc.INFOT = 1;
    zsysv_aa('/', 0, 0, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZSYSV_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zsysv_aa('U', -1, 0, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZSYSV_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zsysv_aa('U', 0, -1, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZSYSV_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zsysv_aa('U', 2, 0, A, 1, IP, B.asMatrix(), 2, W, 1, INFO);
    chkxer('ZSYSV_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zsysv_aa('U', 2, 0, A, 2, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZSYSV_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zsysv_aa('U', 3, 1, A, 3, IP, B.asMatrix(), 3, W, 6, INFO);
    chkxer('ZSYSV_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  } else if (lsamen(2, C2, 'S2')) {
    // ZSYSV_AASEN_2STAGE

    srnamc.SRNAMT = 'ZSYSV_AA_2STAGE';
    infoc.INFOT = 1;
    zsysv_aa_2stage(
        '/', 0, 0, A, 1, A.asArray(), 1, IP, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZSYSV_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zsysv_aa_2stage(
        'U', -1, 0, A, 1, A.asArray(), 1, IP, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZSYSV_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zsysv_aa_2stage(
        'U', 0, -1, A, 1, A.asArray(), 1, IP, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZSYSV_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zsysv_aa_2stage(
        'U', 2, 1, A, 1, A.asArray(), 1, IP, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZSYSV_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zsysv_aa_2stage(
        'U', 2, 1, A, 2, A.asArray(), 1, IP, IP, B.asMatrix(), 2, W, 1, INFO);
    chkxer('ZSYSV_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zsysv_aa_2stage(
        'U', 2, 1, A, 2, A.asArray(), 8, IP, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZSYSV_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    zsysv_aa_2stage(
        'U', 2, 1, A, 2, A.asArray(), 8, IP, IP, B.asMatrix(), 2, W, 1, INFO);
    chkxer('ZSYSV_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  } else if (lsamen(2, C2, 'HP')) {
    // ZHPSV

    srnamc.SRNAMT = 'ZHPSV';
    infoc.INFOT = 1;
    zhpsv('/', 0, 0, A.asArray(), IP, B.asMatrix(), 1, INFO);
    chkxer('ZHPSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zhpsv('U', -1, 0, A.asArray(), IP, B.asMatrix(), 1, INFO);
    chkxer('ZHPSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zhpsv('U', 0, -1, A.asArray(), IP, B.asMatrix(), 1, INFO);
    chkxer('ZHPSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zhpsv('U', 2, 0, A.asArray(), IP, B.asMatrix(), 1, INFO);
    chkxer('ZHPSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZHPSVX

    srnamc.SRNAMT = 'ZHPSVX';
    infoc.INFOT = 1;
    final RCOND = Box(0.0);
    zhpsvx('/', 'U', 0, 0, A.asArray(), AF.asArray(), IP, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZHPSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zhpsvx('N', '/', 0, 0, A.asArray(), AF.asArray(), IP, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZHPSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zhpsvx('N', 'U', -1, 0, A.asArray(), AF.asArray(), IP, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZHPSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zhpsvx('N', 'U', 0, -1, A.asArray(), AF.asArray(), IP, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZHPSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    zhpsvx('N', 'U', 2, 0, A.asArray(), AF.asArray(), IP, B.asMatrix(), 1,
        X.asMatrix(), 2, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZHPSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zhpsvx('N', 'U', 2, 0, A.asArray(), AF.asArray(), IP, B.asMatrix(), 2,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZHPSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  } else if (lsamen(2, C2, 'SY')) {
    // ZSYSV

    srnamc.SRNAMT = 'ZSYSV';
    infoc.INFOT = 1;
    zsysv('/', 0, 0, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZSYSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zsysv('U', -1, 0, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZSYSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zsysv('U', 0, -1, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZSYSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zsysv('U', 2, 0, A, 1, IP, B.asMatrix(), 2, W, 1, INFO);
    chkxer('ZSYSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zsysv('U', 2, 0, A, 2, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZSYSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zsysv('U', 0, 0, A, 1, IP, B.asMatrix(), 1, W, 0, INFO);
    chkxer('ZSYSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zsysv('U', 0, 0, A, 1, IP, B.asMatrix(), 1, W, -2, INFO);
    chkxer('ZSYSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZSYSVX

    srnamc.SRNAMT = 'ZSYSVX';
    infoc.INFOT = 1;
    final RCOND = Box(0.0);
    zsysvx('/', 'U', 0, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
        RCOND, R1, R2, W, 1, RW, INFO);
    chkxer('ZSYSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zsysvx('N', '/', 0, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
        RCOND, R1, R2, W, 1, RW, INFO);
    chkxer('ZSYSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zsysvx('N', 'U', -1, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
        RCOND, R1, R2, W, 1, RW, INFO);
    chkxer('ZSYSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zsysvx('N', 'U', 0, -1, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
        RCOND, R1, R2, W, 1, RW, INFO);
    chkxer('ZSYSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zsysvx('N', 'U', 2, 0, A, 1, AF, 2, IP, B.asMatrix(), 2, X.asMatrix(), 2,
        RCOND, R1, R2, W, 4, RW, INFO);
    chkxer('ZSYSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zsysvx('N', 'U', 2, 0, A, 2, AF, 1, IP, B.asMatrix(), 2, X.asMatrix(), 2,
        RCOND, R1, R2, W, 4, RW, INFO);
    chkxer('ZSYSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zsysvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B.asMatrix(), 1, X.asMatrix(), 2,
        RCOND, R1, R2, W, 4, RW, INFO);
    chkxer('ZSYSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    zsysvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B.asMatrix(), 2, X.asMatrix(), 1,
        RCOND, R1, R2, W, 4, RW, INFO);
    chkxer('ZSYSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 18;
    zsysvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B.asMatrix(), 2, X.asMatrix(), 2,
        RCOND, R1, R2, W, 3, RW, INFO);
    chkxer('ZSYSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  } else if (lsamen(2, C2, 'SR')) {
    // ZSYSV_ROOK

    srnamc.SRNAMT = 'ZSYSV_ROOK';
    infoc.INFOT = 1;
    zsysv_rook('/', 0, 0, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZSYSV_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zsysv_rook('U', -1, 0, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZSYSV_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zsysv_rook('U', 0, -1, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZSYSV_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zsysv_rook('U', 2, 0, A, 1, IP, B.asMatrix(), 2, W, 1, INFO);
    chkxer('ZSYSV_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zsysv_rook('U', 2, 0, A, 2, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZSYSV_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zsysv_rook('U', 0, 0, A, 1, IP, B.asMatrix(), 1, W, 0, INFO);
    chkxer('ZSYSV_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zsysv_rook('U', 0, 0, A, 1, IP, B.asMatrix(), 1, W, -2, INFO);
  } else if (lsamen(2, C2, 'SK')) {
    // ZSYSV_RK

    // Test error exits of the driver that uses factorization
    // of a symmetric indefinite matrix with rook
    // (bounded Bunch-Kaufman) pivoting with the new storage
    // format for factors L ( or U) and D.

    // L (or U) is stored in A, diagonal of D is stored on the
    // diagonal of A, subdiagonal of D is stored in a separate array E.

    srnamc.SRNAMT = 'ZSYSV_RK';
    infoc.INFOT = 1;
    zsysv_rk('/', 0, 0, A, 1, E, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZSYSV_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zsysv_rk('U', -1, 0, A, 1, E, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZSYSV_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zsysv_rk('U', 0, -1, A, 1, E, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZSYSV_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zsysv_rk('U', 2, 0, A, 1, E, IP, B.asMatrix(), 2, W, 1, INFO);
    chkxer('ZSYSV_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    zsysv_rk('U', 2, 0, A, 2, E, IP, B.asMatrix(), 1, W, 1, INFO);
    chkxer('ZSYSV_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zsysv_rk('U', 0, 0, A, 1, E, IP, B.asMatrix(), 1, W, 0, INFO);
    chkxer('ZSYSV_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zsysv_rk('U', 0, 0, A, 1, E, IP, B.asMatrix(), 1, W, -2, INFO);
    chkxer('ZSYSV_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  } else if (lsamen(2, C2, 'SP')) {
    // ZSPSV

    srnamc.SRNAMT = 'ZSPSV';
    infoc.INFOT = 1;
    zspsv('/', 0, 0, A.asArray(), IP, B.asMatrix(), 1, INFO);
    chkxer('ZSPSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zspsv('U', -1, 0, A.asArray(), IP, B.asMatrix(), 1, INFO);
    chkxer('ZSPSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zspsv('U', 0, -1, A.asArray(), IP, B.asMatrix(), 1, INFO);
    chkxer('ZSPSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zspsv('U', 2, 0, A.asArray(), IP, B.asMatrix(), 1, INFO);
    chkxer('ZSPSV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZSPSVX

    srnamc.SRNAMT = 'ZSPSVX';
    infoc.INFOT = 1;
    final RCOND = Box(0.0);
    zspsvx('/', 'U', 0, 0, A.asArray(), AF.asArray(), IP, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZSPSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zspsvx('N', '/', 0, 0, A.asArray(), AF.asArray(), IP, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZSPSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zspsvx('N', 'U', -1, 0, A.asArray(), AF.asArray(), IP, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZSPSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zspsvx('N', 'U', 0, -1, A.asArray(), AF.asArray(), IP, B.asMatrix(), 1,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZSPSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    zspsvx('N', 'U', 2, 0, A.asArray(), AF.asArray(), IP, B.asMatrix(), 1,
        X.asMatrix(), 2, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZSPSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zspsvx('N', 'U', 2, 0, A.asArray(), AF.asArray(), IP, B.asMatrix(), 2,
        X.asMatrix(), 1, RCOND, R1, R2, W, RW, INFO);
    chkxer('ZSPSVX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  }

  // Print a summary line.

  if (infoc.OK.value) {
    NOUT.println(' ${PATH.a3} drivers passed the tests of the error exits');
  } else {
    NOUT.println(
        ' *** ${PATH.a3} drivers failed the tests of the error exits ***');
  }
}
