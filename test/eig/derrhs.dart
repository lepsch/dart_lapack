import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgebak.dart';
import 'package:lapack/src/dgebal.dart';
import 'package:lapack/src/dgehd2.dart';
import 'package:lapack/src/dgehrd.dart';
import 'package:lapack/src/dhsein.dart';
import 'package:lapack/src/dhseqr.dart';
import 'package:lapack/src/dorghr.dart';
import 'package:lapack/src/dormhr.dart';
import 'package:lapack/src/dtrevc.dart';
import 'package:lapack/src/dtrevc3.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../test_driver.dart';
import 'chkxer.dart';
import 'common.dart';

void derrhs(final String PATH, final Nout NUNIT, final TestDriver test) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 3, LW = (NMAX + 2) * (NMAX + 2) + NMAX;
  final SEL = Array<bool>(NMAX);
  final IFAILL = Array<int>(NMAX), IFAILR = Array<int>(NMAX);
  final A = Matrix<double>(NMAX, NMAX),
      C = Matrix<double>(NMAX, NMAX),
      VL = Matrix<double>(NMAX, NMAX),
      VR = Matrix<double>(NMAX, NMAX);
  final S = Array<double>(NMAX),
      TAU = Array<double>(NMAX),
      W = Array<double>(LW),
      WI = Array<double>(NMAX),
      WR = Array<double>(NMAX);
  final INFO = Box(0), M = Box(0), IHI = Box(0), ILO = Box(0);

  final NOUT = infoc.NOUT = NUNIT;
  final OK = infoc.OK;
  final LERR = infoc.LERR;
  NOUT.println();
  final C2 = PATH.substring(1, 3);

  // Set the variables to innocuous values.

  for (var J = 1; J <= NMAX; J++) {
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = 1.0 / (I + J);
    }
    WI[J] = J.toDouble();
    SEL[J] = true;
  }
  OK.value = true;
  var NT = 0;

  // Test error exits of the nonsymmetric eigenvalue routines.

  if (lsamen(2, C2, 'HS')) {
    test('DGEBAL', () {
      srnamc.SRNAMT = 'DGEBAL';
      infoc.INFOT = 1;
      dgebal('/', 0, A, 1, ILO, IHI, S, INFO);
      chkxer('DGEBAL', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgebal('N', -1, A, 1, ILO, IHI, S, INFO);
      chkxer('DGEBAL', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dgebal('N', 2, A, 1, ILO, IHI, S, INFO);
      chkxer('DGEBAL', infoc.INFOT, NOUT, LERR, OK, test);
      NT += 3;
    });

    test('DGEBAK', () {
      srnamc.SRNAMT = 'DGEBAK';
      infoc.INFOT = 1;
      dgebak('/', 'R', 0, 1, 0, S, 0, A, 1, INFO);
      chkxer('DGEBAK', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgebak('N', '/', 0, 1, 0, S, 0, A, 1, INFO);
      chkxer('DGEBAK', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dgebak('N', 'R', -1, 1, 0, S, 0, A, 1, INFO);
      chkxer('DGEBAK', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dgebak('N', 'R', 0, 0, 0, S, 0, A, 1, INFO);
      chkxer('DGEBAK', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dgebak('N', 'R', 0, 2, 0, S, 0, A, 1, INFO);
      chkxer('DGEBAK', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dgebak('N', 'R', 2, 2, 1, S, 0, A, 2, INFO);
      chkxer('DGEBAK', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dgebak('N', 'R', 0, 1, 1, S, 0, A, 1, INFO);
      chkxer('DGEBAK', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 7;
      dgebak('N', 'R', 0, 1, 0, S, -1, A, 1, INFO);
      chkxer('DGEBAK', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 9;
      dgebak('N', 'R', 2, 1, 2, S, 0, A, 1, INFO);
      chkxer('DGEBAK', infoc.INFOT, NOUT, LERR, OK, test);
      NT += 9;
    });

    test('DGEHRD', () {
      srnamc.SRNAMT = 'DGEHRD';
      infoc.INFOT = 1;
      dgehrd(-1, 1, 1, A, 1, TAU, W, 1, INFO);
      chkxer('DGEHRD', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgehrd(0, 0, 0, A, 1, TAU, W, 1, INFO);
      chkxer('DGEHRD', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgehrd(0, 2, 0, A, 1, TAU, W, 1, INFO);
      chkxer('DGEHRD', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dgehrd(1, 1, 0, A, 1, TAU, W, 1, INFO);
      chkxer('DGEHRD', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dgehrd(0, 1, 1, A, 1, TAU, W, 1, INFO);
      chkxer('DGEHRD', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dgehrd(2, 1, 1, A, 1, TAU, W, 2, INFO);
      chkxer('DGEHRD', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      dgehrd(2, 1, 2, A, 2, TAU, W, 1, INFO);
      chkxer('DGEHRD', infoc.INFOT, NOUT, LERR, OK, test);
      NT += 7;
    });

    test('DGEHD2', () {
      srnamc.SRNAMT = 'DGEHD2';
      infoc.INFOT = 1;
      dgehd2(-1, 1, 1, A, 1, TAU, W, INFO);
      chkxer('DGEHD2', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgehd2(0, 0, 0, A, 1, TAU, W, INFO);
      chkxer('DGEHD2', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dgehd2(0, 2, 0, A, 1, TAU, W, INFO);
      chkxer('DGEHD2', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dgehd2(1, 1, 0, A, 1, TAU, W, INFO);
      chkxer('DGEHD2', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dgehd2(0, 1, 1, A, 1, TAU, W, INFO);
      chkxer('DGEHD2', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dgehd2(2, 1, 1, A, 1, TAU, W, INFO);
      chkxer('DGEHD2', infoc.INFOT, NOUT, LERR, OK, test);
      NT += 6;
    });

    test('DORGHR', () {
      srnamc.SRNAMT = 'DORGHR';
      infoc.INFOT = 1;
      dorghr(-1, 1, 1, A, 1, TAU, W, 1, INFO);
      chkxer('DORGHR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dorghr(0, 0, 0, A, 1, TAU, W, 1, INFO);
      chkxer('DORGHR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dorghr(0, 2, 0, A, 1, TAU, W, 1, INFO);
      chkxer('DORGHR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dorghr(1, 1, 0, A, 1, TAU, W, 1, INFO);
      chkxer('DORGHR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dorghr(0, 1, 1, A, 1, TAU, W, 1, INFO);
      chkxer('DORGHR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dorghr(2, 1, 1, A, 1, TAU, W, 1, INFO);
      chkxer('DORGHR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      dorghr(3, 1, 3, A, 3, TAU, W, 1, INFO);
      chkxer('DORGHR', infoc.INFOT, NOUT, LERR, OK, test);
      NT += 7;
    });

    test('DORMHR', () {
      srnamc.SRNAMT = 'DORMHR';
      infoc.INFOT = 1;
      dormhr('/', 'N', 0, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO);
      chkxer('DORMHR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dormhr('L', '/', 0, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO);
      chkxer('DORMHR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dormhr('L', 'N', -1, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO);
      chkxer('DORMHR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dormhr('L', 'N', 0, -1, 1, 0, A, 1, TAU, C, 1, W, 1, INFO);
      chkxer('DORMHR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dormhr('L', 'N', 0, 0, 0, 0, A, 1, TAU, C, 1, W, 1, INFO);
      chkxer('DORMHR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dormhr('L', 'N', 0, 0, 2, 0, A, 1, TAU, C, 1, W, 1, INFO);
      chkxer('DORMHR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dormhr('L', 'N', 1, 2, 2, 1, A, 1, TAU, C, 1, W, 2, INFO);
      chkxer('DORMHR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dormhr('R', 'N', 2, 1, 2, 1, A, 1, TAU, C, 2, W, 2, INFO);
      chkxer('DORMHR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dormhr('L', 'N', 1, 1, 1, 0, A, 1, TAU, C, 1, W, 1, INFO);
      chkxer('DORMHR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dormhr('L', 'N', 0, 1, 1, 1, A, 1, TAU, C, 1, W, 1, INFO);
      chkxer('DORMHR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dormhr('R', 'N', 1, 0, 1, 1, A, 1, TAU, C, 1, W, 1, INFO);
      chkxer('DORMHR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      dormhr('L', 'N', 2, 1, 1, 1, A, 1, TAU, C, 2, W, 1, INFO);
      chkxer('DORMHR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      dormhr('R', 'N', 1, 2, 1, 1, A, 1, TAU, C, 1, W, 1, INFO);
      chkxer('DORMHR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 11;
      dormhr('L', 'N', 2, 1, 1, 1, A, 2, TAU, C, 1, W, 1, INFO);
      chkxer('DORMHR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 13;
      dormhr('L', 'N', 1, 2, 1, 1, A, 1, TAU, C, 1, W, 1, INFO);
      chkxer('DORMHR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 13;
      dormhr('R', 'N', 2, 1, 1, 1, A, 1, TAU, C, 2, W, 1, INFO);
      chkxer('DORMHR', infoc.INFOT, NOUT, LERR, OK, test);
      NT += 16;
    });

    test('DHSEQR', () {
      srnamc.SRNAMT = 'DHSEQR';
      infoc.INFOT = 1;
      dhseqr('/', 'N', 0, 1, 0, A, 1, WR, WI, C, 1, W, 1, INFO);
      chkxer('DHSEQR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dhseqr('E', '/', 0, 1, 0, A, 1, WR, WI, C, 1, W, 1, INFO);
      chkxer('DHSEQR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dhseqr('E', 'N', -1, 1, 0, A, 1, WR, WI, C, 1, W, 1, INFO);
      chkxer('DHSEQR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dhseqr('E', 'N', 0, 0, 0, A, 1, WR, WI, C, 1, W, 1, INFO);
      chkxer('DHSEQR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dhseqr('E', 'N', 0, 2, 0, A, 1, WR, WI, C, 1, W, 1, INFO);
      chkxer('DHSEQR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dhseqr('E', 'N', 1, 1, 0, A, 1, WR, WI, C, 1, W, 1, INFO);
      chkxer('DHSEQR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dhseqr('E', 'N', 1, 1, 2, A, 1, WR, WI, C, 1, W, 1, INFO);
      chkxer('DHSEQR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 7;
      dhseqr('E', 'N', 2, 1, 2, A, 1, WR, WI, C, 2, W, 1, INFO);
      chkxer('DHSEQR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 11;
      dhseqr('E', 'V', 2, 1, 2, A, 2, WR, WI, C, 1, W, 1, INFO);
      chkxer('DHSEQR', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 13;
      dhseqr('E', 'N', 2, 1, 2, A, 2, WR, WI, C, 1, W, 1, INFO);
      chkxer('DHSEQR', infoc.INFOT, NOUT, LERR, OK, test);
      NT += 10;
    });

    test('DHSEIN', () {
      srnamc.SRNAMT = 'DHSEIN';
      infoc.INFOT = 1;
      dhsein('/', 'N', 'N', SEL, 0, A, 1, WR, WI, VL, 1, VR, 1, 0, M, W, IFAILL,
          IFAILR, INFO);
      chkxer('DHSEIN', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dhsein('R', '/', 'N', SEL, 0, A, 1, WR, WI, VL, 1, VR, 1, 0, M, W, IFAILL,
          IFAILR, INFO);
      chkxer('DHSEIN', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dhsein('R', 'N', '/', SEL, 0, A, 1, WR, WI, VL, 1, VR, 1, 0, M, W, IFAILL,
          IFAILR, INFO);
      chkxer('DHSEIN', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dhsein('R', 'N', 'N', SEL, -1, A, 1, WR, WI, VL, 1, VR, 1, 0, M, W,
          IFAILL, IFAILR, INFO);
      chkxer('DHSEIN', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 7;
      dhsein('R', 'N', 'N', SEL, 2, A, 1, WR, WI, VL, 1, VR, 2, 4, M, W, IFAILL,
          IFAILR, INFO);
      chkxer('DHSEIN', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 11;
      dhsein('L', 'N', 'N', SEL, 2, A, 2, WR, WI, VL, 1, VR, 1, 4, M, W, IFAILL,
          IFAILR, INFO);
      chkxer('DHSEIN', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 13;
      dhsein('R', 'N', 'N', SEL, 2, A, 2, WR, WI, VL, 1, VR, 1, 4, M, W, IFAILL,
          IFAILR, INFO);
      chkxer('DHSEIN', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 14;
      dhsein('R', 'N', 'N', SEL, 2, A, 2, WR, WI, VL, 1, VR, 2, 1, M, W, IFAILL,
          IFAILR, INFO);
      chkxer('DHSEIN', infoc.INFOT, NOUT, LERR, OK, test);
      NT += 8;
    });

    test('DTREVC', () {
      srnamc.SRNAMT = 'DTREVC';
      infoc.INFOT = 1;
      dtrevc('/', 'A', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, INFO);
      chkxer('DTREVC', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dtrevc('L', '/', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, INFO);
      chkxer('DTREVC', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dtrevc('L', 'A', SEL, -1, A, 1, VL, 1, VR, 1, 0, M, W, INFO);
      chkxer('DTREVC', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dtrevc('L', 'A', SEL, 2, A, 1, VL, 2, VR, 1, 4, M, W, INFO);
      chkxer('DTREVC', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      dtrevc('L', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, INFO);
      chkxer('DTREVC', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 10;
      dtrevc('R', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, INFO);
      chkxer('DTREVC', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 11;
      dtrevc('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 1, M, W, INFO);
      chkxer('DTREVC', infoc.INFOT, NOUT, LERR, OK, test);
      NT += 7;
    });

    test('DTREVC3', () {
      srnamc.SRNAMT = 'DTREVC3';
      infoc.INFOT = 1;
      dtrevc3('/', 'A', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, LW, INFO);
      chkxer('DTREVC3', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dtrevc3('L', '/', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, LW, INFO);
      chkxer('DTREVC3', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dtrevc3('L', 'A', SEL, -1, A, 1, VL, 1, VR, 1, 0, M, W, LW, INFO);
      chkxer('DTREVC3', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dtrevc3('L', 'A', SEL, 2, A, 1, VL, 2, VR, 1, 4, M, W, LW, INFO);
      chkxer('DTREVC3', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      dtrevc3('L', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, LW, INFO);
      chkxer('DTREVC3', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 10;
      dtrevc3('R', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, LW, INFO);
      chkxer('DTREVC3', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 11;
      dtrevc3('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 1, M, W, LW, INFO);
      chkxer('DTREVC3', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 14;
      dtrevc3('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 2, M, W, 2, INFO);
      chkxer('DTREVC3', infoc.INFOT, NOUT, LERR, OK, test);
      NT += 8;
    });
  }

  // Print a summary line.

  if (OK.value) {
    NOUT.println(
        ' ${PATH.a3} routines passed the tests of the error exits (${NT.i3} tests done)');
  } else {
    NOUT.println(
        ' *** ${PATH.a3} routines failed the tests of the error exits ***');
  }
}
