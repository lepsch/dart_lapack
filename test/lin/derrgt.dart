import 'package:lapack/lapack.dart';

import '../test_driver.dart';
import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void derrgt(final String PATH, final Nout NUNIT, final TestDriver test) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 2;
  final IP = Array<int>(NMAX), IW = Array<int>(NMAX);
  final B = Array<double>(NMAX),
      C = Array<double>(NMAX),
      CF = Array<double>(NMAX),
      D = Array<double>(NMAX),
      DF = Array<double>(NMAX),
      E = Array<double>(NMAX),
      EF = Array<double>(NMAX),
      F = Array<double>(NMAX),
      R1 = Array<double>(NMAX),
      R2 = Array<double>(NMAX),
      W = Array<double>(NMAX),
      X = Array<double>(NMAX);
  final INFO = Box(0);
  final ANORM = Box(0.0), RCOND = Box(0.0);

  final NOUT = infoc.NOUT = NUNIT;
  final OK = infoc.OK;
  final LERR = infoc.LERR;
  NOUT.println();
  final C2 = PATH.substring(1, 3);
  D[1] = 1.0;
  D[2] = 2.0;
  DF[1] = 1.0;
  DF[2] = 2.0;
  E[1] = 3.0;
  E[2] = 4.0;
  EF[1] = 3.0;
  EF[2] = 4.0;
  ANORM.value = 1.0;
  OK.value = true;

  if (lsamen(2, C2, 'GT')) {
    // Test error exits for the general tridiagonal routines.

    test('DGTTRF', () {
      srnamc.SRNAMT = 'DGTTRF';
      infoc.INFOT = 1;
      dgttrf(-1, C, D, E, F, IP, INFO);
      chkxer('DGTTRF', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DGTTRS', () {
      srnamc.SRNAMT = 'DGTTRS';
      infoc.INFOT = 1;
      dgttrs('/', 0, 0, C, D, E, F, IP, X.asMatrix(), 1, INFO);
      chkxer('DGTTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgttrs('N', -1, 0, C, D, E, F, IP, X.asMatrix(), 1, INFO);
      chkxer('DGTTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dgttrs('N', 0, -1, C, D, E, F, IP, X.asMatrix(), 1, INFO);
      chkxer('DGTTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 10;
      dgttrs('N', 2, 1, C, D, E, F, IP, X.asMatrix(), 1, INFO);
      chkxer('DGTTRS', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DGTRFS', () {
      srnamc.SRNAMT = 'DGTRFS';
      infoc.INFOT = 1;
      dgtrfs('/', 0, 0, C, D, E, CF, DF, EF, F, IP, B.asMatrix(), 1,
          X.asMatrix(), 1, R1, R2, W, IW, INFO);
      chkxer('DGTRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgtrfs('N', -1, 0, C, D, E, CF, DF, EF, F, IP, B.asMatrix(), 1,
          X.asMatrix(), 1, R1, R2, W, IW, INFO);
      chkxer('DGTRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dgtrfs('N', 0, -1, C, D, E, CF, DF, EF, F, IP, B.asMatrix(), 1,
          X.asMatrix(), 1, R1, R2, W, IW, INFO);
      chkxer('DGTRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 13;
      dgtrfs('N', 2, 1, C, D, E, CF, DF, EF, F, IP, B.asMatrix(), 1,
          X.asMatrix(), 2, R1, R2, W, IW, INFO);
      chkxer('DGTRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 15;
      dgtrfs('N', 2, 1, C, D, E, CF, DF, EF, F, IP, B.asMatrix(), 2,
          X.asMatrix(), 1, R1, R2, W, IW, INFO);
      chkxer('DGTRFS', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DGTCON', () {
      srnamc.SRNAMT = 'DGTCON';
      infoc.INFOT = 1;
      dgtcon('/', 0, C, D, E, F, IP, ANORM.value, RCOND, W, IW, INFO);
      chkxer('DGTCON', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgtcon('I', -1, C, D, E, F, IP, ANORM.value, RCOND, W, IW, INFO);
      chkxer('DGTCON', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      dgtcon('I', 0, C, D, E, F, IP, -ANORM.value, RCOND, W, IW, INFO);
      chkxer('DGTCON', infoc.INFOT, NOUT, LERR, OK, test);
    });
  } else if (lsamen(2, C2, 'PT')) {
    // Test error exits for the positive definite tridiagonal
    // routines.

    test('DPTTRF', () {
      srnamc.SRNAMT = 'DPTTRF';
      infoc.INFOT = 1;
      dpttrf(-1, D, E, INFO);
      chkxer('DPTTRF', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DPTTRS', () {
      srnamc.SRNAMT = 'DPTTRS';
      infoc.INFOT = 1;
      dpttrs(-1, 0, D, E, X.asMatrix(), 1, INFO);
      chkxer('DPTTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dpttrs(0, -1, D, E, X.asMatrix(), 1, INFO);
      chkxer('DPTTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dpttrs(2, 1, D, E, X.asMatrix(), 1, INFO);
      chkxer('DPTTRS', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DPTRFS', () {
      srnamc.SRNAMT = 'DPTRFS';
      infoc.INFOT = 1;
      dptrfs(-1, 0, D, E, DF, EF, B.asMatrix(), 1, X.asMatrix(), 1, R1, R2, W,
          INFO);
      chkxer('DPTRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dptrfs(0, -1, D, E, DF, EF, B.asMatrix(), 1, X.asMatrix(), 1, R1, R2, W,
          INFO);
      chkxer('DPTRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      dptrfs(2, 1, D, E, DF, EF, B.asMatrix(), 1, X.asMatrix(), 2, R1, R2, W,
          INFO);
      chkxer('DPTRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 10;
      dptrfs(2, 1, D, E, DF, EF, B.asMatrix(), 2, X.asMatrix(), 1, R1, R2, W,
          INFO);
      chkxer('DPTRFS', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DPTCON', () {
      srnamc.SRNAMT = 'DPTCON';
      infoc.INFOT = 1;
      dptcon(-1, D, E, ANORM.value, RCOND, W, INFO);
      chkxer('DPTCON', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dptcon(0, D, E, -ANORM.value, RCOND, W, INFO);
      chkxer('DPTCON', infoc.INFOT, NOUT, LERR, OK, test);
    });
  }

  // Print a summary line.
  alaesm(PATH, OK.value, NOUT);
}
