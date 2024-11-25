// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/lsamen.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/zgbcon.dart';
import 'package:dart_lapack/src/zgbequ.dart';
import 'package:dart_lapack/src/zgbequb.dart';
import 'package:dart_lapack/src/zgbrfs.dart';
import 'package:dart_lapack/src/zgbrfsx.dart';
import 'package:dart_lapack/src/zgbtf2.dart';
import 'package:dart_lapack/src/zgbtrf.dart';
import 'package:dart_lapack/src/zgbtrs.dart';
import 'package:dart_lapack/src/zgecon.dart';
import 'package:dart_lapack/src/zgeequ.dart';
import 'package:dart_lapack/src/zgeequb.dart';
import 'package:dart_lapack/src/zgerfs.dart';
import 'package:dart_lapack/src/zgerfsx.dart';
import 'package:dart_lapack/src/zgetf2.dart';
import 'package:dart_lapack/src/zgetrf.dart';
import 'package:dart_lapack/src/zgetri.dart';
import 'package:dart_lapack/src/zgetrs.dart';

import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void zerrge(final String PATH, final Nout NUNIT) {
  const NMAX = 4;
  String EQ;
  int N_ERR_BNDS, NPARAMS;
  final IP = Array<int>(NMAX);
  final R = Array<double>(NMAX),
      R1 = Array<double>(NMAX),
      R2 = Array<double>(NMAX),
      CS = Array<double>(NMAX),
      RS = Array<double>(NMAX),
      BERR = Array<double>(1),
      PARAMS = Array<double>(1);
  final A = Matrix<Complex>(NMAX, NMAX),
      AF = Matrix<Complex>(NMAX, NMAX),
      B = Array<Complex>(NMAX),
      W = Array<Complex>(2 * NMAX),
      X = Array<Complex>(NMAX),
      ERR_BNDS_N = Matrix<Complex>(NMAX, 3),
      ERR_BNDS_C = Matrix<Complex>(NMAX, 3);
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
    R1[J] = 0.0;
    R2[J] = 0.0;
    W[J] = Complex.zero;
    X[J] = Complex.zero;
    CS[J] = 0.0;
    RS[J] = 0.0;
    IP[J] = J;
  }
  infoc.OK.value = true;

  // Test error exits of the routines that use the LU decomposition
  // of a general matrix.

  if (lsamen(2, C2, 'GE')) {
    // ZGETRF

    srnamc.SRNAMT = 'ZGETRF';
    infoc.INFOT = 1;
    zgetrf(-1, 0, A, 1, IP, INFO);
    chkxer('ZGETRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgetrf(0, -1, A, 1, IP, INFO);
    chkxer('ZGETRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgetrf(2, 1, A, 1, IP, INFO);
    chkxer('ZGETRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGETF2

    srnamc.SRNAMT = 'ZGETF2';
    infoc.INFOT = 1;
    zgetf2(-1, 0, A, 1, IP, INFO);
    chkxer('ZGETF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgetf2(0, -1, A, 1, IP, INFO);
    chkxer('ZGETF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgetf2(2, 1, A, 1, IP, INFO);
    chkxer('ZGETF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGETRI

    srnamc.SRNAMT = 'ZGETRI';
    infoc.INFOT = 1;
    zgetri(-1, A, 1, IP, W, 1, INFO);
    chkxer('ZGETRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgetri(2, A, 1, IP, W, 2, INFO);
    chkxer('ZGETRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zgetri(2, A, 2, IP, W, 1, INFO);
    chkxer('ZGETRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGETRS

    srnamc.SRNAMT = 'ZGETRS';
    infoc.INFOT = 1;
    zgetrs('/', 0, 0, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('ZGETRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgetrs('N', -1, 0, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('ZGETRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgetrs('N', 0, -1, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('ZGETRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgetrs('N', 2, 1, A, 1, IP, B.asMatrix(), 2, INFO);
    chkxer('ZGETRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zgetrs('N', 2, 1, A, 2, IP, B.asMatrix(), 1, INFO);
    chkxer('ZGETRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGERFS

    srnamc.SRNAMT = 'ZGERFS';
    infoc.INFOT = 1;
    zgerfs('/', 0, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1, R1, R2,
        W, R, INFO);
    chkxer('ZGERFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgerfs('N', -1, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1, R1,
        R2, W, R, INFO);
    chkxer('ZGERFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgerfs('N', 0, -1, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1, R1,
        R2, W, R, INFO);
    chkxer('ZGERFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgerfs('N', 2, 1, A, 1, AF, 2, IP, B.asMatrix(), 2, X.asMatrix(), 2, R1, R2,
        W, R, INFO);
    chkxer('ZGERFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zgerfs('N', 2, 1, A, 2, AF, 1, IP, B.asMatrix(), 2, X.asMatrix(), 2, R1, R2,
        W, R, INFO);
    chkxer('ZGERFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zgerfs('N', 2, 1, A, 2, AF, 2, IP, B.asMatrix(), 1, X.asMatrix(), 2, R1, R2,
        W, R, INFO);
    chkxer('ZGERFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    zgerfs('N', 2, 1, A, 2, AF, 2, IP, B.asMatrix(), 2, X.asMatrix(), 1, R1, R2,
        W, R, INFO);
    chkxer('ZGERFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGERFSX

    N_ERR_BNDS = 3;
    NPARAMS = 0;
    srnamc.SRNAMT = 'ZGERFSX';
    infoc.INFOT = 1;
    EQ = ' ';
    final RCOND = Box(0.0);
    zgerfsx(
        '/',
        EQ,
        0,
        0,
        A,
        1,
        AF,
        1,
        IP,
        RS,
        CS,
        B.asMatrix(),
        1,
        X.asMatrix(),
        1,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N.cast<double>(),
        ERR_BNDS_C.cast<double>(),
        NPARAMS,
        PARAMS,
        W,
        R,
        INFO);
    chkxer('ZGERFSX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    EQ = '/';
    zgerfsx(
        'N',
        EQ,
        2,
        1,
        A,
        1,
        AF,
        2,
        IP,
        RS,
        CS,
        B.asMatrix(),
        2,
        X.asMatrix(),
        2,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N.cast<double>(),
        ERR_BNDS_C.cast<double>(),
        NPARAMS,
        PARAMS,
        W,
        R,
        INFO);
    chkxer('ZGERFSX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    EQ = 'R';
    zgerfsx(
        'N',
        EQ,
        -1,
        0,
        A,
        1,
        AF,
        1,
        IP,
        RS,
        CS,
        B.asMatrix(),
        1,
        X.asMatrix(),
        1,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N.cast<double>(),
        ERR_BNDS_C.cast<double>(),
        NPARAMS,
        PARAMS,
        W,
        R,
        INFO);
    chkxer('ZGERFSX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgerfsx(
        'N',
        EQ,
        0,
        -1,
        A,
        1,
        AF,
        1,
        IP,
        RS,
        CS,
        B.asMatrix(),
        1,
        X.asMatrix(),
        1,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N.cast<double>(),
        ERR_BNDS_C.cast<double>(),
        NPARAMS,
        PARAMS,
        W,
        R,
        INFO);
    chkxer('ZGERFSX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zgerfsx(
        'N',
        EQ,
        2,
        1,
        A,
        1,
        AF,
        2,
        IP,
        RS,
        CS,
        B.asMatrix(),
        2,
        X.asMatrix(),
        2,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N.cast<double>(),
        ERR_BNDS_C.cast<double>(),
        NPARAMS,
        PARAMS,
        W,
        R,
        INFO);
    chkxer('ZGERFSX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zgerfsx(
        'N',
        EQ,
        2,
        1,
        A,
        2,
        AF,
        1,
        IP,
        RS,
        CS,
        B.asMatrix(),
        2,
        X.asMatrix(),
        2,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N.cast<double>(),
        ERR_BNDS_C.cast<double>(),
        NPARAMS,
        PARAMS,
        W,
        R,
        INFO);
    chkxer('ZGERFSX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    EQ = 'C';
    zgerfsx(
        'N',
        EQ,
        2,
        1,
        A,
        2,
        AF,
        2,
        IP,
        RS,
        CS,
        B.asMatrix(),
        1,
        X.asMatrix(),
        2,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N.cast<double>(),
        ERR_BNDS_C.cast<double>(),
        NPARAMS,
        PARAMS,
        W,
        R,
        INFO);
    chkxer('ZGERFSX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 15;
    zgerfsx(
        'N',
        EQ,
        2,
        1,
        A,
        2,
        AF,
        2,
        IP,
        RS,
        CS,
        B.asMatrix(),
        2,
        X.asMatrix(),
        1,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N.cast<double>(),
        ERR_BNDS_C.cast<double>(),
        NPARAMS,
        PARAMS,
        W,
        R,
        INFO);
    chkxer('ZGERFSX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGECON

    srnamc.SRNAMT = 'ZGECON';
    infoc.INFOT = 1;
    final ANRM = Box(0.0);
    zgecon('/', 0, A, 1, ANRM.value, RCOND, W, R, INFO);
    chkxer('ZGECON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgecon('1', -1, A, 1, ANRM.value, RCOND, W, R, INFO);
    chkxer('ZGECON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgecon('1', 2, A, 1, ANRM.value, RCOND, W, R, INFO);
    chkxer('ZGECON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGEEQU

    srnamc.SRNAMT = 'ZGEEQU';
    infoc.INFOT = 1;
    final CCOND = Box(0.0);
    zgeequ(-1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('ZGEEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgeequ(0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('ZGEEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgeequ(2, 2, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('ZGEEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGEEQUB

    srnamc.SRNAMT = 'ZGEEQUB';
    infoc.INFOT = 1;
    zgeequb(-1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('ZGEEQUB', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgeequb(0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('ZGEEQUB', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgeequb(2, 2, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('ZGEEQUB', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // Test error exits of the routines that use the LU decomposition
    // of a general band matrix.
  } else if (lsamen(2, C2, 'GB')) {
    // ZGBTRF

    srnamc.SRNAMT = 'ZGBTRF';
    infoc.INFOT = 1;
    zgbtrf(-1, 0, 0, 0, A, 1, IP, INFO);
    chkxer('ZGBTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgbtrf(0, -1, 0, 0, A, 1, IP, INFO);
    chkxer('ZGBTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgbtrf(1, 1, -1, 0, A, 1, IP, INFO);
    chkxer('ZGBTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgbtrf(1, 1, 0, -1, A, 1, IP, INFO);
    chkxer('ZGBTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zgbtrf(2, 2, 1, 1, A, 3, IP, INFO);
    chkxer('ZGBTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGBTF2

    srnamc.SRNAMT = 'ZGBTF2';
    infoc.INFOT = 1;
    zgbtf2(-1, 0, 0, 0, A, 1, IP, INFO);
    chkxer('ZGBTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgbtf2(0, -1, 0, 0, A, 1, IP, INFO);
    chkxer('ZGBTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgbtf2(1, 1, -1, 0, A, 1, IP, INFO);
    chkxer('ZGBTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgbtf2(1, 1, 0, -1, A, 1, IP, INFO);
    chkxer('ZGBTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zgbtf2(2, 2, 1, 1, A, 3, IP, INFO);
    chkxer('ZGBTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGBTRS

    srnamc.SRNAMT = 'ZGBTRS';
    infoc.INFOT = 1;
    zgbtrs('/', 0, 0, 0, 1, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('ZGBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgbtrs('N', -1, 0, 0, 1, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('ZGBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgbtrs('N', 1, -1, 0, 1, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('ZGBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgbtrs('N', 1, 0, -1, 1, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('ZGBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgbtrs('N', 1, 0, 0, -1, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('ZGBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zgbtrs('N', 2, 1, 1, 1, A, 3, IP, B.asMatrix(), 2, INFO);
    chkxer('ZGBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zgbtrs('N', 2, 0, 0, 1, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('ZGBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGBRFS

    srnamc.SRNAMT = 'ZGBRFS';
    infoc.INFOT = 1;
    zgbrfs('/', 0, 0, 0, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
        R1, R2, W, R, INFO);
    chkxer('ZGBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgbrfs('N', -1, 0, 0, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
        R1, R2, W, R, INFO);
    chkxer('ZGBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgbrfs('N', 1, -1, 0, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
        R1, R2, W, R, INFO);
    chkxer('ZGBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgbrfs('N', 1, 0, -1, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
        R1, R2, W, R, INFO);
    chkxer('ZGBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgbrfs('N', 1, 0, 0, -1, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
        R1, R2, W, R, INFO);
    chkxer('ZGBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zgbrfs('N', 2, 1, 1, 1, A, 2, AF, 4, IP, B.asMatrix(), 2, X.asMatrix(), 2,
        R1, R2, W, R, INFO);
    chkxer('ZGBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    zgbrfs('N', 2, 1, 1, 1, A, 3, AF, 3, IP, B.asMatrix(), 2, X.asMatrix(), 2,
        R1, R2, W, R, INFO);
    chkxer('ZGBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    zgbrfs('N', 2, 0, 0, 1, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 2,
        R1, R2, W, R, INFO);
    chkxer('ZGBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 14;
    zgbrfs('N', 2, 0, 0, 1, A, 1, AF, 1, IP, B.asMatrix(), 2, X.asMatrix(), 1,
        R1, R2, W, R, INFO);
    chkxer('ZGBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGBRFSX

    N_ERR_BNDS = 3;
    NPARAMS = 0;
    srnamc.SRNAMT = 'ZGBRFSX';
    infoc.INFOT = 1;
    EQ = ' ';
    final RCOND = Box(0.0);
    zgbrfsx(
        '/',
        EQ,
        0,
        0,
        0,
        0,
        A,
        1,
        AF,
        1,
        IP,
        RS,
        CS,
        B.asMatrix(),
        1,
        X.asMatrix(),
        1,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N.cast<double>(),
        ERR_BNDS_C.cast<double>(),
        NPARAMS,
        PARAMS,
        W,
        R,
        INFO);
    chkxer('ZGBRFSX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    EQ = '/';
    zgbrfsx(
        'N',
        EQ,
        2,
        1,
        1,
        1,
        A,
        1,
        AF,
        2,
        IP,
        RS,
        CS,
        B.asMatrix(),
        2,
        X.asMatrix(),
        2,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N.cast<double>(),
        ERR_BNDS_C.cast<double>(),
        NPARAMS,
        PARAMS,
        W,
        R,
        INFO);
    chkxer('ZGBRFSX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    EQ = 'R';
    zgbrfsx(
        'N',
        EQ,
        -1,
        1,
        1,
        0,
        A,
        1,
        AF,
        1,
        IP,
        RS,
        CS,
        B.asMatrix(),
        1,
        X.asMatrix(),
        1,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N.cast<double>(),
        ERR_BNDS_C.cast<double>(),
        NPARAMS,
        PARAMS,
        W,
        R,
        INFO);
    chkxer('ZGBRFSX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    EQ = 'R';
    zgbrfsx(
        'N',
        EQ,
        2,
        -1,
        1,
        1,
        A,
        3,
        AF,
        4,
        IP,
        RS,
        CS,
        B.asMatrix(),
        1,
        X.asMatrix(),
        1,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N.cast<double>(),
        ERR_BNDS_C.cast<double>(),
        NPARAMS,
        PARAMS,
        W,
        R,
        INFO);
    chkxer('ZGBRFSX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    EQ = 'R';
    zgbrfsx(
        'N',
        EQ,
        2,
        1,
        -1,
        1,
        A,
        3,
        AF,
        4,
        IP,
        RS,
        CS,
        B.asMatrix(),
        1,
        X.asMatrix(),
        1,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N.cast<double>(),
        ERR_BNDS_C.cast<double>(),
        NPARAMS,
        PARAMS,
        W,
        R,
        INFO);
    chkxer('ZGBRFSX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zgbrfsx(
        'N',
        EQ,
        0,
        0,
        0,
        -1,
        A,
        1,
        AF,
        1,
        IP,
        RS,
        CS,
        B.asMatrix(),
        1,
        X.asMatrix(),
        1,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N.cast<double>(),
        ERR_BNDS_C.cast<double>(),
        NPARAMS,
        PARAMS,
        W,
        R,
        INFO);
    chkxer('ZGBRFSX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zgbrfsx(
        'N',
        EQ,
        2,
        1,
        1,
        1,
        A,
        1,
        AF,
        2,
        IP,
        RS,
        CS,
        B.asMatrix(),
        2,
        X.asMatrix(),
        2,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N.cast<double>(),
        ERR_BNDS_C.cast<double>(),
        NPARAMS,
        PARAMS,
        W,
        R,
        INFO);
    chkxer('ZGBRFSX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zgbrfsx(
        'N',
        EQ,
        2,
        1,
        1,
        1,
        A,
        3,
        AF,
        3,
        IP,
        RS,
        CS,
        B.asMatrix(),
        2,
        X.asMatrix(),
        2,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N.cast<double>(),
        ERR_BNDS_C.cast<double>(),
        NPARAMS,
        PARAMS,
        W,
        R,
        INFO);
    chkxer('ZGBRFSX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    EQ = 'C';
    zgbrfsx(
        'N',
        EQ,
        2,
        1,
        1,
        1,
        A,
        3,
        AF,
        5,
        IP,
        RS,
        CS,
        B.asMatrix(),
        1,
        X.asMatrix(),
        2,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N.cast<double>(),
        ERR_BNDS_C.cast<double>(),
        NPARAMS,
        PARAMS,
        W,
        R,
        INFO);
    chkxer('ZGBRFSX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 15;
    zgbrfsx(
        'N',
        EQ,
        2,
        1,
        1,
        1,
        A,
        3,
        AF,
        5,
        IP,
        RS,
        CS,
        B.asMatrix(),
        2,
        X.asMatrix(),
        1,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N.cast<double>(),
        ERR_BNDS_C.cast<double>(),
        NPARAMS,
        PARAMS,
        W,
        R,
        INFO);
    chkxer('ZGBRFSX', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGBCON

    srnamc.SRNAMT = 'ZGBCON';
    infoc.INFOT = 1;
    final ANRM = Box(0.0);
    zgbcon('/', 0, 0, 0, A, 1, IP, ANRM.value, RCOND, W, R, INFO);
    chkxer('ZGBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgbcon('1', -1, 0, 0, A, 1, IP, ANRM.value, RCOND, W, R, INFO);
    chkxer('ZGBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgbcon('1', 1, -1, 0, A, 1, IP, ANRM.value, RCOND, W, R, INFO);
    chkxer('ZGBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgbcon('1', 1, 0, -1, A, 1, IP, ANRM.value, RCOND, W, R, INFO);
    chkxer('ZGBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zgbcon('1', 2, 1, 1, A, 3, IP, ANRM.value, RCOND, W, R, INFO);
    chkxer('ZGBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGBEQU

    srnamc.SRNAMT = 'ZGBEQU';
    infoc.INFOT = 1;
    final CCOND = Box(0.0);
    zgbequ(-1, 0, 0, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('ZGBEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgbequ(0, -1, 0, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('ZGBEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgbequ(1, 1, -1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('ZGBEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgbequ(1, 1, 0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('ZGBEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zgbequ(2, 2, 1, 1, A, 2, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('ZGBEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZGBEQUB

    srnamc.SRNAMT = 'ZGBEQUB';
    infoc.INFOT = 1;
    zgbequb(-1, 0, 0, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('ZGBEQUB', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgbequb(0, -1, 0, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('ZGBEQUB', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgbequb(1, 1, -1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('ZGBEQUB', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgbequb(1, 1, 0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('ZGBEQUB', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zgbequb(2, 2, 1, 1, A, 2, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('ZGBEQUB', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  }

  // Print a summary line.

  alaesm(PATH, infoc.OK.value, NOUT);
}
