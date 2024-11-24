// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zgbcon.dart';
import 'package:lapack/src/zgbequ.dart';
import 'package:lapack/src/zgbrfs.dart';
import 'package:lapack/src/zgbtf2.dart';
import 'package:lapack/src/zgbtrf.dart';
import 'package:lapack/src/zgbtrs.dart';
import 'package:lapack/src/zgecon.dart';
import 'package:lapack/src/zgeequ.dart';
import 'package:lapack/src/zgerfs.dart';
import 'package:lapack/src/zgetf2.dart';
import 'package:lapack/src/zgetrf.dart';
import 'package:lapack/src/zgetri.dart';
import 'package:lapack/src/zgetrs.dart';

import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void zerrge(final String PATH, final Nout NUNIT) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 4;
  final IP = Array<int>(NMAX);
  final R = Array<double>(NMAX),
      R1 = Array<double>(NMAX),
      R2 = Array<double>(NMAX);
  final A = Matrix<Complex>(NMAX, NMAX),
      AF = Matrix<Complex>(NMAX, NMAX),
      B = Array<Complex>(NMAX),
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
    R1[J] = 0.0;
    R2[J] = 0.0;
    W[J] = Complex.zero;
    X[J] = Complex.zero;
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

    // ZGECON

    srnamc.SRNAMT = 'ZGECON';
    infoc.INFOT = 1;
    final RCOND = Box(0.0), ANRM = Box(0.0), CCOND = Box(0.0);
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
    zgeequ(-1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('ZGEEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgeequ(0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('ZGEEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgeequ(2, 2, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('ZGEEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

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

    // ZGBCON

    srnamc.SRNAMT = 'ZGBCON';
    infoc.INFOT = 1;
    final RCOND = Box(0.0), ANRM = Box(0.0), CCOND = Box(0.0);
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
  }

  // Print a summary line.

  alaesm(PATH, infoc.OK.value, NOUT);
}
