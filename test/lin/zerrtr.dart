import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zlatbs.dart';
import 'package:lapack/src/zlatps.dart';
import 'package:lapack/src/zlatrs.dart';
import 'package:lapack/src/zlatrs3.dart';
import 'package:lapack/src/ztbcon.dart';
import 'package:lapack/src/ztbrfs.dart';
import 'package:lapack/src/ztbtrs.dart';
import 'package:lapack/src/ztpcon.dart';
import 'package:lapack/src/ztprfs.dart';
import 'package:lapack/src/ztptri.dart';
import 'package:lapack/src/ztptrs.dart';
import 'package:lapack/src/ztrcon.dart';
import 'package:lapack/src/ztrrfs.dart';
import 'package:lapack/src/ztrti2.dart';
import 'package:lapack/src/ztrtri.dart';
import 'package:lapack/src/ztrtrs.dart';

import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void zerrtr(final String PATH, final Nout NUNIT) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 2;
  final SCALES = Array<double>(0),
      R1 = Array<double>(NMAX),
      R2 = Array<double>(NMAX),
      RW = Array<double>(NMAX);
  final A = Matrix<Complex>(NMAX, NMAX),
      B = Array<Complex>(NMAX),
      W = Array<Complex>(NMAX),
      X = Array<Complex>(NMAX);
  final INFO = Box(0);

  final NOUT = NUNIT;
  NOUT.println();
  final C2 = PATH.substring(1, 3);
  A[1][1] = 1.0.toComplex();
  A[1][2] = 2.0.toComplex();
  A[2][2] = 3.0.toComplex();
  A[2][1] = 4.0.toComplex();
  infoc.OK.value = true;

  // Test error exits for the general triangular routines.

  if (lsamen(2, C2, 'TR')) {
    // ZTRTRI

    srnamc.SRNAMT = 'ZTRTRI';
    infoc.INFOT = 1;
    ztrtri('/', 'N', 0, A, 1, INFO);
    chkxer('ZTRTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    ztrtri('U', '/', 0, A, 1, INFO);
    chkxer('ZTRTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    ztrtri('U', 'N', -1, A, 1, INFO);
    chkxer('ZTRTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    ztrtri('U', 'N', 2, A, 1, INFO);
    chkxer('ZTRTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZTRTI2

    srnamc.SRNAMT = 'ZTRTI2';
    infoc.INFOT = 1;
    ztrti2('/', 'N', 0, A, 1, INFO);
    chkxer('ZTRTI2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    ztrti2('U', '/', 0, A, 1, INFO);
    chkxer('ZTRTI2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    ztrti2('U', 'N', -1, A, 1, INFO);
    chkxer('ZTRTI2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    ztrti2('U', 'N', 2, A, 1, INFO);
    chkxer('ZTRTI2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZTRTRS

    srnamc.SRNAMT = 'ZTRTRS';
    infoc.INFOT = 1;
    ztrtrs('/', 'N', 'N', 0, 0, A, 1, X.asMatrix(), 1, INFO);
    chkxer('ZTRTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    ztrtrs('U', '/', 'N', 0, 0, A, 1, X.asMatrix(), 1, INFO);
    chkxer('ZTRTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    ztrtrs('U', 'N', '/', 0, 0, A, 1, X.asMatrix(), 1, INFO);
    chkxer('ZTRTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    ztrtrs('U', 'N', 'N', -1, 0, A, 1, X.asMatrix(), 1, INFO);
    chkxer('ZTRTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    ztrtrs('U', 'N', 'N', 0, -1, A, 1, X.asMatrix(), 1, INFO);
    chkxer('ZTRTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;

    // ZTRRFS

    srnamc.SRNAMT = 'ZTRRFS';
    infoc.INFOT = 1;
    ztrrfs('/', 'N', 'N', 0, 0, A, 1, B.asMatrix(), 1, X.asMatrix(), 1, R1, R2,
        W, RW, INFO);
    chkxer('ZTRRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    ztrrfs('U', '/', 'N', 0, 0, A, 1, B.asMatrix(), 1, X.asMatrix(), 1, R1, R2,
        W, RW, INFO);
    chkxer('ZTRRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    ztrrfs('U', 'N', '/', 0, 0, A, 1, B.asMatrix(), 1, X.asMatrix(), 1, R1, R2,
        W, RW, INFO);
    chkxer('ZTRRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    ztrrfs('U', 'N', 'N', -1, 0, A, 1, B.asMatrix(), 1, X.asMatrix(), 1, R1, R2,
        W, RW, INFO);
    chkxer('ZTRRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    ztrrfs('U', 'N', 'N', 0, -1, A, 1, B.asMatrix(), 1, X.asMatrix(), 1, R1, R2,
        W, RW, INFO);
    chkxer('ZTRRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    ztrrfs('U', 'N', 'N', 2, 1, A, 1, B.asMatrix(), 2, X.asMatrix(), 2, R1, R2,
        W, RW, INFO);
    chkxer('ZTRRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    ztrrfs('U', 'N', 'N', 2, 1, A, 2, B.asMatrix(), 1, X.asMatrix(), 2, R1, R2,
        W, RW, INFO);
    chkxer('ZTRRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    ztrrfs('U', 'N', 'N', 2, 1, A, 2, B.asMatrix(), 2, X.asMatrix(), 1, R1, R2,
        W, RW, INFO);
    chkxer('ZTRRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZTRCON

    srnamc.SRNAMT = 'ZTRCON';
    infoc.INFOT = 1;
    final RCOND = Box(0.0);
    ztrcon('/', 'U', 'N', 0, A, 1, RCOND, W, RW, INFO);
    chkxer('ZTRCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    ztrcon('1', '/', 'N', 0, A, 1, RCOND, W, RW, INFO);
    chkxer('ZTRCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    ztrcon('1', 'U', '/', 0, A, 1, RCOND, W, RW, INFO);
    chkxer('ZTRCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    ztrcon('1', 'U', 'N', -1, A, 1, RCOND, W, RW, INFO);
    chkxer('ZTRCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    ztrcon('1', 'U', 'N', 2, A, 1, RCOND, W, RW, INFO);
    chkxer('ZTRCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZLATRS

    srnamc.SRNAMT = 'ZLATRS';
    infoc.INFOT = 1;
    final SCALE = Box(0.0);
    zlatrs('/', 'N', 'N', 'N', 0, A, 1, X, SCALE, RW, INFO);
    chkxer('ZLATRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zlatrs('U', '/', 'N', 'N', 0, A, 1, X, SCALE, RW, INFO);
    chkxer('ZLATRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zlatrs('U', 'N', '/', 'N', 0, A, 1, X, SCALE, RW, INFO);
    chkxer('ZLATRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zlatrs('U', 'N', 'N', '/', 0, A, 1, X, SCALE, RW, INFO);
    chkxer('ZLATRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zlatrs('U', 'N', 'N', 'N', -1, A, 1, X, SCALE, RW, INFO);
    chkxer('ZLATRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zlatrs('U', 'N', 'N', 'N', 2, A, 1, X, SCALE, RW, INFO);
    chkxer('ZLATRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZLATRS3

    srnamc.SRNAMT = 'ZLATRS3';
    infoc.INFOT = 1;
    zlatrs3('/', 'N', 'N', 'N', 0, 0, A, 1, X.asMatrix(), 1, SCALES, RW, RW(2),
        1, INFO);
    chkxer('ZLATRS3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zlatrs3('U', '/', 'N', 'N', 0, 0, A, 1, X.asMatrix(), 1, SCALES, RW, RW(2),
        1, INFO);
    chkxer('ZLATRS3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zlatrs3('U', 'N', '/', 'N', 0, 0, A, 1, X.asMatrix(), 1, SCALES, RW, RW(2),
        1, INFO);
    chkxer('ZLATRS3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zlatrs3('U', 'N', 'N', '/', 0, 0, A, 1, X.asMatrix(), 1, SCALES, RW, RW(2),
        1, INFO);
    chkxer('ZLATRS3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zlatrs3('U', 'N', 'N', 'N', -1, 0, A, 1, X.asMatrix(), 1, SCALES, RW, RW(2),
        1, INFO);
    chkxer('ZLATRS3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zlatrs3('U', 'N', 'N', 'N', 0, -1, A, 1, X.asMatrix(), 1, SCALES, RW, RW(2),
        1, INFO);
    chkxer('ZLATRS3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zlatrs3('U', 'N', 'N', 'N', 2, 0, A, 1, X.asMatrix(), 1, SCALES, RW, RW(2),
        1, INFO);
    chkxer('ZLATRS3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zlatrs3('U', 'N', 'N', 'N', 2, 0, A, 2, X.asMatrix(), 1, SCALES, RW, RW(2),
        1, INFO);
    chkxer('ZLATRS3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 14;
    zlatrs3('U', 'N', 'N', 'N', 1, 0, A, 1, X.asMatrix(), 1, SCALES, RW, RW(2),
        0, INFO);
    chkxer('ZLATRS3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // Test error exits for the packed triangular routines.
  } else if (lsamen(2, C2, 'TP')) {
    // ZTPTRI

    srnamc.SRNAMT = 'ZTPTRI';
    infoc.INFOT = 1;
    ztptri('/', 'N', 0, A.asArray(), INFO);
    chkxer('ZTPTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    ztptri('U', '/', 0, A.asArray(), INFO);
    chkxer('ZTPTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    ztptri('U', 'N', -1, A.asArray(), INFO);
    chkxer('ZTPTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZTPTRS

    srnamc.SRNAMT = 'ZTPTRS';
    infoc.INFOT = 1;
    ztptrs('/', 'N', 'N', 0, 0, A.asArray(), X.asMatrix(), 1, INFO);
    chkxer('ZTPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    ztptrs('U', '/', 'N', 0, 0, A.asArray(), X.asMatrix(), 1, INFO);
    chkxer('ZTPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    ztptrs('U', 'N', '/', 0, 0, A.asArray(), X.asMatrix(), 1, INFO);
    chkxer('ZTPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    ztptrs('U', 'N', 'N', -1, 0, A.asArray(), X.asMatrix(), 1, INFO);
    chkxer('ZTPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    ztptrs('U', 'N', 'N', 0, -1, A.asArray(), X.asMatrix(), 1, INFO);
    chkxer('ZTPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    ztptrs('U', 'N', 'N', 2, 1, A.asArray(), X.asMatrix(), 1, INFO);
    chkxer('ZTPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZTPRFS

    srnamc.SRNAMT = 'ZTPRFS';
    infoc.INFOT = 1;
    ztprfs('/', 'N', 'N', 0, 0, A.asArray(), B.asMatrix(), 1, X.asMatrix(), 1,
        R1, R2, W, RW, INFO);
    chkxer('ZTPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    ztprfs('U', '/', 'N', 0, 0, A.asArray(), B.asMatrix(), 1, X.asMatrix(), 1,
        R1, R2, W, RW, INFO);
    chkxer('ZTPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    ztprfs('U', 'N', '/', 0, 0, A.asArray(), B.asMatrix(), 1, X.asMatrix(), 1,
        R1, R2, W, RW, INFO);
    chkxer('ZTPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    ztprfs('U', 'N', 'N', -1, 0, A.asArray(), B.asMatrix(), 1, X.asMatrix(), 1,
        R1, R2, W, RW, INFO);
    chkxer('ZTPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    ztprfs('U', 'N', 'N', 0, -1, A.asArray(), B.asMatrix(), 1, X.asMatrix(), 1,
        R1, R2, W, RW, INFO);
    chkxer('ZTPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    ztprfs('U', 'N', 'N', 2, 1, A.asArray(), B.asMatrix(), 1, X.asMatrix(), 2,
        R1, R2, W, RW, INFO);
    chkxer('ZTPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    ztprfs('U', 'N', 'N', 2, 1, A.asArray(), B.asMatrix(), 2, X.asMatrix(), 1,
        R1, R2, W, RW, INFO);
    chkxer('ZTPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZTPCON

    srnamc.SRNAMT = 'ZTPCON';
    infoc.INFOT = 1;
    final RCOND = Box(0.0);
    ztpcon('/', 'U', 'N', 0, A.asArray(), RCOND, W, RW, INFO);
    chkxer('ZTPCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    ztpcon('1', '/', 'N', 0, A.asArray(), RCOND, W, RW, INFO);
    chkxer('ZTPCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    ztpcon('1', 'U', '/', 0, A.asArray(), RCOND, W, RW, INFO);
    chkxer('ZTPCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    ztpcon('1', 'U', 'N', -1, A.asArray(), RCOND, W, RW, INFO);
    chkxer('ZTPCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZLATPS

    srnamc.SRNAMT = 'ZLATPS';
    infoc.INFOT = 1;
    final SCALE = Box(0.0);
    zlatps('/', 'N', 'N', 'N', 0, A.asArray(), X, SCALE, RW, INFO);
    chkxer('ZLATPS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zlatps('U', '/', 'N', 'N', 0, A.asArray(), X, SCALE, RW, INFO);
    chkxer('ZLATPS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zlatps('U', 'N', '/', 'N', 0, A.asArray(), X, SCALE, RW, INFO);
    chkxer('ZLATPS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zlatps('U', 'N', 'N', '/', 0, A.asArray(), X, SCALE, RW, INFO);
    chkxer('ZLATPS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zlatps('U', 'N', 'N', 'N', -1, A.asArray(), X, SCALE, RW, INFO);
    chkxer('ZLATPS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // Test error exits for the banded triangular routines.
  } else if (lsamen(2, C2, 'TB')) {
    // ZTBTRS

    srnamc.SRNAMT = 'ZTBTRS';
    infoc.INFOT = 1;
    ztbtrs('/', 'N', 'N', 0, 0, 0, A, 1, X.asMatrix(), 1, INFO);
    chkxer('ZTBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    ztbtrs('U', '/', 'N', 0, 0, 0, A, 1, X.asMatrix(), 1, INFO);
    chkxer('ZTBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    ztbtrs('U', 'N', '/', 0, 0, 0, A, 1, X.asMatrix(), 1, INFO);
    chkxer('ZTBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    ztbtrs('U', 'N', 'N', -1, 0, 0, A, 1, X.asMatrix(), 1, INFO);
    chkxer('ZTBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    ztbtrs('U', 'N', 'N', 0, -1, 0, A, 1, X.asMatrix(), 1, INFO);
    chkxer('ZTBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    ztbtrs('U', 'N', 'N', 0, 0, -1, A, 1, X.asMatrix(), 1, INFO);
    chkxer('ZTBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    ztbtrs('U', 'N', 'N', 2, 1, 1, A, 1, X.asMatrix(), 2, INFO);
    chkxer('ZTBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    ztbtrs('U', 'N', 'N', 2, 0, 1, A, 1, X.asMatrix(), 1, INFO);
    chkxer('ZTBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZTBRFS

    srnamc.SRNAMT = 'ZTBRFS';
    infoc.INFOT = 1;
    ztbrfs('/', 'N', 'N', 0, 0, 0, A, 1, B.asMatrix(), 1, X.asMatrix(), 1, R1,
        R2, W, RW, INFO);
    chkxer('ZTBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    ztbrfs('U', '/', 'N', 0, 0, 0, A, 1, B.asMatrix(), 1, X.asMatrix(), 1, R1,
        R2, W, RW, INFO);
    chkxer('ZTBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    ztbrfs('U', 'N', '/', 0, 0, 0, A, 1, B.asMatrix(), 1, X.asMatrix(), 1, R1,
        R2, W, RW, INFO);
    chkxer('ZTBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    ztbrfs('U', 'N', 'N', -1, 0, 0, A, 1, B.asMatrix(), 1, X.asMatrix(), 1, R1,
        R2, W, RW, INFO);
    chkxer('ZTBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    ztbrfs('U', 'N', 'N', 0, -1, 0, A, 1, B.asMatrix(), 1, X.asMatrix(), 1, R1,
        R2, W, RW, INFO);
    chkxer('ZTBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    ztbrfs('U', 'N', 'N', 0, 0, -1, A, 1, B.asMatrix(), 1, X.asMatrix(), 1, R1,
        R2, W, RW, INFO);
    chkxer('ZTBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    ztbrfs('U', 'N', 'N', 2, 1, 1, A, 1, B.asMatrix(), 2, X.asMatrix(), 2, R1,
        R2, W, RW, INFO);
    chkxer('ZTBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    ztbrfs('U', 'N', 'N', 2, 1, 1, A, 2, B.asMatrix(), 1, X.asMatrix(), 2, R1,
        R2, W, RW, INFO);
    chkxer('ZTBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    ztbrfs('U', 'N', 'N', 2, 1, 1, A, 2, B.asMatrix(), 2, X.asMatrix(), 1, R1,
        R2, W, RW, INFO);
    chkxer('ZTBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZTBCON

    srnamc.SRNAMT = 'ZTBCON';
    infoc.INFOT = 1;
    final RCOND = Box(0.0);
    ztbcon('/', 'U', 'N', 0, 0, A, 1, RCOND, W, RW, INFO);
    chkxer('ZTBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    ztbcon('1', '/', 'N', 0, 0, A, 1, RCOND, W, RW, INFO);
    chkxer('ZTBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    ztbcon('1', 'U', '/', 0, 0, A, 1, RCOND, W, RW, INFO);
    chkxer('ZTBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    ztbcon('1', 'U', 'N', -1, 0, A, 1, RCOND, W, RW, INFO);
    chkxer('ZTBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    ztbcon('1', 'U', 'N', 0, -1, A, 1, RCOND, W, RW, INFO);
    chkxer('ZTBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    ztbcon('1', 'U', 'N', 2, 1, A, 1, RCOND, W, RW, INFO);
    chkxer('ZTBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

    // ZLATBS

    srnamc.SRNAMT = 'ZLATBS';
    infoc.INFOT = 1;
    final SCALE = Box(0.0);
    zlatbs('/', 'N', 'N', 'N', 0, 0, A, 1, X, SCALE, RW, INFO);
    chkxer('ZLATBS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zlatbs('U', '/', 'N', 'N', 0, 0, A, 1, X, SCALE, RW, INFO);
    chkxer('ZLATBS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zlatbs('U', 'N', '/', 'N', 0, 0, A, 1, X, SCALE, RW, INFO);
    chkxer('ZLATBS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zlatbs('U', 'N', 'N', '/', 0, 0, A, 1, X, SCALE, RW, INFO);
    chkxer('ZLATBS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zlatbs('U', 'N', 'N', 'N', -1, 0, A, 1, X, SCALE, RW, INFO);
    chkxer('ZLATBS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zlatbs('U', 'N', 'N', 'N', 1, -1, A, 1, X, SCALE, RW, INFO);
    chkxer('ZLATBS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zlatbs('U', 'N', 'N', 'N', 2, 1, A, 1, X, SCALE, RW, INFO);
    chkxer('ZLATBS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  }

  // Print a summary line.

  alaesm(PATH, infoc.OK.value, NOUT);
}
