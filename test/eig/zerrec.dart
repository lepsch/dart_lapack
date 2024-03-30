import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/ztrexc.dart';
import 'package:lapack/src/ztrsen.dart';
import 'package:lapack/src/ztrsna.dart';
import 'package:lapack/src/ztrsyl.dart';
import 'package:lapack/src/ztrsyl3.dart';

import 'chkxer.dart';
import 'common.dart';

void zerrec(
  final String PATH,
  final Nout NUNIT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final NMAX = Box(4);
  final LW = NMAX.value * (NMAX.value + 2);
  int I, IFST, ILST, J, NT;
  final SCALE = Box(0.0);
  final INFO = Box(0);
  final SEL = Array<bool>(NMAX.value);
  final RW = Array<double>(LW),
      S = Array<double>(NMAX.value),
      SEP = Array<double>(NMAX.value),
      SWORK = Array<double>(NMAX.value);
  final A = Matrix<Complex>(NMAX.value, NMAX.value),
      B = Matrix<Complex>(NMAX.value, NMAX.value),
      C = Matrix<Complex>(NMAX.value, NMAX.value),
      WORK = Array<Complex>(LW),
      X = Array<Complex>(NMAX.value);
  final M = Box(0);

  infoc.NOUT = NUNIT;
  infoc.OK.value = true;
  NT = 0;

  // Initialize A, B and SEL

  for (J = 1; J <= NMAX.value; J++) {
    for (I = 1; I <= NMAX.value; I++) {
      A[I][J] = Complex.zero;
      B[I][J] = Complex.zero;
    }
  }
  for (I = 1; I <= NMAX.value; I++) {
    A[I][I] = Complex.one;
    SEL[I] = true;
  }

  // Test ZTRSYL

  srnamc.SRNAMT = 'ZTRSYL';
  infoc.INFOT = 1;
  ztrsyl('X', 'N', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO);
  chkxer('ZTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  ztrsyl('N', 'X', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO);
  chkxer('ZTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  ztrsyl('N', 'N', 0, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO);
  chkxer('ZTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  ztrsyl('N', 'N', 1, -1, 0, A, 1, B, 1, C, 1, SCALE, INFO);
  chkxer('ZTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  ztrsyl('N', 'N', 1, 0, -1, A, 1, B, 1, C, 1, SCALE, INFO);
  chkxer('ZTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  ztrsyl('N', 'N', 1, 2, 0, A, 1, B, 1, C, 2, SCALE, INFO);
  chkxer('ZTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 9;
  ztrsyl('N', 'N', 1, 0, 2, A, 1, B, 1, C, 1, SCALE, INFO);
  chkxer('ZTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 11;
  ztrsyl('N', 'N', 1, 2, 0, A, 2, B, 1, C, 1, SCALE, INFO);
  chkxer('ZTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  NT += 8;

  // Test ZTRSYL3

  srnamc.SRNAMT = 'ZTRSYL3';
  infoc.INFOT = 1;
  ztrsyl3('X', 'N', 1, 0, 0, A, 1, B, 1, C, 1, SCALE,
      SWORK.asMatrix(NMAX.value), NMAX, INFO);
  chkxer('ZTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  ztrsyl3('N', 'X', 1, 0, 0, A, 1, B, 1, C, 1, SCALE,
      SWORK.asMatrix(NMAX.value), NMAX, INFO);
  chkxer('ZTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  ztrsyl3('N', 'N', 0, 0, 0, A, 1, B, 1, C, 1, SCALE,
      SWORK.asMatrix(NMAX.value), NMAX, INFO);
  chkxer('ZTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  ztrsyl3('N', 'N', 1, -1, 0, A, 1, B, 1, C, 1, SCALE,
      SWORK.asMatrix(NMAX.value), NMAX, INFO);
  chkxer('ZTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  ztrsyl3('N', 'N', 1, 0, -1, A, 1, B, 1, C, 1, SCALE,
      SWORK.asMatrix(NMAX.value), NMAX, INFO);
  chkxer('ZTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  ztrsyl3('N', 'N', 1, 2, 0, A, 1, B, 1, C, 2, SCALE,
      SWORK.asMatrix(NMAX.value), NMAX, INFO);
  chkxer('ZTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 9;
  ztrsyl3('N', 'N', 1, 0, 2, A, 1, B, 1, C, 1, SCALE,
      SWORK.asMatrix(NMAX.value), NMAX, INFO);
  chkxer('ZTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 11;
  ztrsyl3('N', 'N', 1, 2, 0, A, 2, B, 1, C, 1, SCALE,
      SWORK.asMatrix(NMAX.value), NMAX, INFO);
  chkxer('ZTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  NT += 8;

  // Test ZTREXC

  srnamc.SRNAMT = 'ZTREXC';
  IFST = 1;
  ILST = 1;
  infoc.INFOT = 1;
  ztrexc('X', 1, A, 1, B, 1, IFST, ILST, INFO);
  chkxer('ZTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  ztrexc('N', -1, A, 1, B, 1, IFST, ILST, INFO);
  chkxer('ZTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  ILST = 2;
  ztrexc('N', 2, A, 1, B, 1, IFST, ILST, INFO);
  chkxer('ZTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 6;
  ztrexc('V', 2, A, 2, B, 1, IFST, ILST, INFO);
  chkxer('ZTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  IFST = 0;
  ILST = 1;
  ztrexc('V', 1, A, 1, B, 1, IFST, ILST, INFO);
  chkxer('ZTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  IFST = 2;
  ztrexc('V', 1, A, 1, B, 1, IFST, ILST, INFO);
  chkxer('ZTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  IFST = 1;
  ILST = 0;
  ztrexc('V', 1, A, 1, B, 1, IFST, ILST, INFO);
  chkxer('ZTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  ILST = 2;
  ztrexc('V', 1, A, 1, B, 1, IFST, ILST, INFO);
  chkxer('ZTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  NT += 8;

  // Test ZTRSNA

  srnamc.SRNAMT = 'ZTRSNA';
  infoc.INFOT = 1;
  ztrsna('X', 'A', SEL, 0, A, 1, B, 1, C, 1, S, SEP, 1, M, WORK.asMatrix(), 1,
      RW, INFO);
  chkxer('ZTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  ztrsna('B', 'X', SEL, 0, A, 1, B, 1, C, 1, S, SEP, 1, M, WORK.asMatrix(), 1,
      RW, INFO);
  chkxer('ZTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  ztrsna('B', 'A', SEL, -1, A, 1, B, 1, C, 1, S, SEP, 1, M, WORK.asMatrix(), 1,
      RW, INFO);
  chkxer('ZTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 6;
  ztrsna('V', 'A', SEL, 2, A, 1, B, 1, C, 1, S, SEP, 2, M, WORK.asMatrix(), 2,
      RW, INFO);
  chkxer('ZTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  ztrsna('B', 'A', SEL, 2, A, 2, B, 1, C, 2, S, SEP, 2, M, WORK.asMatrix(), 2,
      RW, INFO);
  chkxer('ZTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 10;
  ztrsna('B', 'A', SEL, 2, A, 2, B, 2, C, 1, S, SEP, 2, M, WORK.asMatrix(), 2,
      RW, INFO);
  chkxer('ZTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 13;
  ztrsna('B', 'A', SEL, 1, A, 1, B, 1, C, 1, S, SEP, 0, M, WORK.asMatrix(), 1,
      RW, INFO);
  chkxer('ZTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 13;
  ztrsna('B', 'S', SEL, 2, A, 2, B, 2, C, 2, S, SEP, 1, M, WORK.asMatrix(), 1,
      RW, INFO);
  chkxer('ZTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 16;
  ztrsna('B', 'A', SEL, 2, A, 2, B, 2, C, 2, S, SEP, 2, M, WORK.asMatrix(), 1,
      RW, INFO);
  chkxer('ZTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  NT += 9;

  // Test ZTRSEN

  SEL[1] = false;
  srnamc.SRNAMT = 'ZTRSEN';
  infoc.INFOT = 1;
  ztrsen('X', 'N', SEL, 0, A, 1, B, 1, X, M, S(1), SEP(1), WORK, 1, INFO);
  chkxer('ZTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  ztrsen('N', 'X', SEL, 0, A, 1, B, 1, X, M, S(1), SEP(1), WORK, 1, INFO);
  chkxer('ZTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  ztrsen('N', 'N', SEL, -1, A, 1, B, 1, X, M, S(1), SEP(1), WORK, 1, INFO);
  chkxer('ZTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 6;
  ztrsen('N', 'N', SEL, 2, A, 1, B, 1, X, M, S(1), SEP(1), WORK, 2, INFO);
  chkxer('ZTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  ztrsen('N', 'V', SEL, 2, A, 2, B, 1, X, M, S(1), SEP(1), WORK, 1, INFO);
  chkxer('ZTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 14;
  ztrsen('N', 'V', SEL, 2, A, 2, B, 2, X, M, S(1), SEP(1), WORK, 0, INFO);
  chkxer('ZTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 14;
  ztrsen('E', 'V', SEL, 3, A, 3, B, 3, X, M, S(1), SEP(1), WORK, 1, INFO);
  chkxer('ZTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 14;
  ztrsen('V', 'V', SEL, 3, A, 3, B, 3, X, M, S(1), SEP(1), WORK, 3, INFO);
  chkxer('ZTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  NT += 8;

  // Print a summary line.

  if (infoc.OK.value) {
    infoc.NOUT.println(
        ' ${PATH.a3} routines passed the tests of the error exits (${NT.i3} tests done)');
  } else {
    infoc.NOUT.println(
        ' *** ${PATH.a3} routines failed the tests of the error exits ***');
  }
}
