// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dtrexc.dart';
import 'package:dart_lapack/src/dtrsen.dart';
import 'package:dart_lapack/src/dtrsna.dart';
import 'package:dart_lapack/src/dtrsyl.dart';
import 'package:dart_lapack/src/dtrsyl3.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';

import '../test_driver.dart';
import 'chkxer.dart';
import 'common.dart';

void derrec(final String PATH, final Nout NUNIT, final TestDriver test) {
  const NMAX = 4, ONE = 1.0, ZERO = 0.0;
  final SCALE = Box(0.0);
  final SEL = Array<bool>(NMAX);
  final IWORK = Array<int>(NMAX);
  final A = Matrix<double>(NMAX, NMAX),
      B = Matrix<double>(NMAX, NMAX),
      C = Matrix<double>(NMAX, NMAX);
  final S = Array<double>(NMAX),
      SEP = Array<double>(NMAX),
      WI = Array<double>(NMAX),
      WORK = Array<double>(NMAX),
      WR = Array<double>(NMAX);
  final INFO = Box(0), M = Box(0), IFST = Box(0), ILST = Box(0);

  final NOUT = infoc.NOUT = NUNIT;
  final OK = infoc.OK;
  final LERR = infoc.LERR;
  OK.value = true;
  var NT = 0;

  // Initialize A, B and SEL

  for (var J = 1; J <= NMAX; J++) {
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = ZERO;
      B[I][J] = ZERO;
    }
  }
  for (var I = 1; I <= NMAX; I++) {
    A[I][I] = ONE;
    SEL[I] = true;
  }

  test('DTRSYL', () {
    srnamc.SRNAMT = 'DTRSYL';
    infoc.INFOT = 1;
    dtrsyl('X', 'N', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO);
    chkxer('DTRSYL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dtrsyl('N', 'X', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO);
    chkxer('DTRSYL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 3;
    dtrsyl('N', 'N', 0, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO);
    chkxer('DTRSYL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 4;
    dtrsyl('N', 'N', 1, -1, 0, A, 1, B, 1, C, 1, SCALE, INFO);
    chkxer('DTRSYL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    dtrsyl('N', 'N', 1, 0, -1, A, 1, B, 1, C, 1, SCALE, INFO);
    chkxer('DTRSYL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 7;
    dtrsyl('N', 'N', 1, 2, 0, A, 1, B, 1, C, 2, SCALE, INFO);
    chkxer('DTRSYL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 9;
    dtrsyl('N', 'N', 1, 0, 2, A, 1, B, 1, C, 1, SCALE, INFO);
    chkxer('DTRSYL', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 11;
    dtrsyl('N', 'N', 1, 2, 0, A, 2, B, 1, C, 1, SCALE, INFO);
    chkxer('DTRSYL', infoc.INFOT, NOUT, LERR, OK);
    NT += 8;
  });

  test('DTRSYL3', () {
    srnamc.SRNAMT = 'DTRSYL3';
    infoc.INFOT = 1;
    var NMAXINOUT = Box(NMAX);
    dtrsyl3('X', 'N', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX,
        WORK.asMatrix(NMAX), NMAXINOUT, INFO);
    chkxer('DTRSYL3', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    NMAXINOUT = Box(NMAX);
    dtrsyl3('N', 'X', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX,
        WORK.asMatrix(NMAX), NMAXINOUT, INFO);
    chkxer('DTRSYL3', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 3;
    NMAXINOUT = Box(NMAX);
    dtrsyl3('N', 'N', 0, 0, 0, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX,
        WORK.asMatrix(NMAX), NMAXINOUT, INFO);
    chkxer('DTRSYL3', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 4;
    NMAXINOUT = Box(NMAX);
    dtrsyl3('N', 'N', 1, -1, 0, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX,
        WORK.asMatrix(NMAX), NMAXINOUT, INFO);
    chkxer('DTRSYL3', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 5;
    NMAXINOUT = Box(NMAX);
    dtrsyl3('N', 'N', 1, 0, -1, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX,
        WORK.asMatrix(NMAX), NMAXINOUT, INFO);
    chkxer('DTRSYL3', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 7;
    NMAXINOUT = Box(NMAX);
    dtrsyl3('N', 'N', 1, 2, 0, A, 1, B, 1, C, 2, SCALE, IWORK, NMAX,
        WORK.asMatrix(NMAX), NMAXINOUT, INFO);
    chkxer('DTRSYL3', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 9;
    NMAXINOUT = Box(NMAX);
    dtrsyl3('N', 'N', 1, 0, 2, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX,
        WORK.asMatrix(NMAX), NMAXINOUT, INFO);
    chkxer('DTRSYL3', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 11;
    NMAXINOUT = Box(NMAX);
    dtrsyl3('N', 'N', 1, 2, 0, A, 2, B, 1, C, 1, SCALE, IWORK, NMAX,
        WORK.asMatrix(NMAX), NMAXINOUT, INFO);
    chkxer('DTRSYL3', infoc.INFOT, NOUT, LERR, OK);
    NT += 8;
  });

  test('DTREXC', () {
    srnamc.SRNAMT = 'DTREXC';
    IFST.value = 1;
    ILST.value = 1;
    infoc.INFOT = 1;
    dtrexc('X', 1, A, 1, B, 1, IFST, ILST, WORK, INFO);
    chkxer('DTREXC', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dtrexc('N', -1, A, 1, B, 1, IFST, ILST, WORK, INFO);
    chkxer('DTREXC', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 4;
    ILST.value = 2;
    dtrexc('N', 2, A, 1, B, 1, IFST, ILST, WORK, INFO);
    chkxer('DTREXC', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 6;
    dtrexc('V', 2, A, 2, B, 1, IFST, ILST, WORK, INFO);
    chkxer('DTREXC', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 7;
    IFST.value = 0;
    ILST.value = 1;
    dtrexc('V', 1, A, 1, B, 1, IFST, ILST, WORK, INFO);
    chkxer('DTREXC', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 7;
    IFST.value = 2;
    dtrexc('V', 1, A, 1, B, 1, IFST, ILST, WORK, INFO);
    chkxer('DTREXC', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 8;
    IFST.value = 1;
    ILST.value = 0;
    dtrexc('V', 1, A, 1, B, 1, IFST, ILST, WORK, INFO);
    chkxer('DTREXC', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 8;
    ILST.value = 2;
    dtrexc('V', 1, A, 1, B, 1, IFST, ILST, WORK, INFO);
    chkxer('DTREXC', infoc.INFOT, NOUT, LERR, OK);
    NT += 8;
  });

  test('DTRSNA', () {
    srnamc.SRNAMT = 'DTRSNA';
    infoc.INFOT = 1;
    dtrsna('X', 'A', SEL, 0, A, 1, B, 1, C, 1, S, SEP, 1, M, WORK.asMatrix(1),
        1, IWORK, INFO);
    chkxer('DTRSNA', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dtrsna('B', 'X', SEL, 0, A, 1, B, 1, C, 1, S, SEP, 1, M, WORK.asMatrix(1),
        1, IWORK, INFO);
    chkxer('DTRSNA', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 4;
    dtrsna('B', 'A', SEL, -1, A, 1, B, 1, C, 1, S, SEP, 1, M, WORK.asMatrix(1),
        1, IWORK, INFO);
    chkxer('DTRSNA', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 6;
    dtrsna('V', 'A', SEL, 2, A, 1, B, 1, C, 1, S, SEP, 2, M, WORK.asMatrix(2),
        2, IWORK, INFO);
    chkxer('DTRSNA', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 8;
    dtrsna('B', 'A', SEL, 2, A, 2, B, 1, C, 2, S, SEP, 2, M, WORK.asMatrix(2),
        2, IWORK, INFO);
    chkxer('DTRSNA', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 10;
    dtrsna('B', 'A', SEL, 2, A, 2, B, 2, C, 1, S, SEP, 2, M, WORK.asMatrix(2),
        2, IWORK, INFO);
    chkxer('DTRSNA', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 13;
    dtrsna('B', 'A', SEL, 1, A, 1, B, 1, C, 1, S, SEP, 0, M, WORK.asMatrix(1),
        1, IWORK, INFO);
    chkxer('DTRSNA', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 13;
    dtrsna('B', 'S', SEL, 2, A, 2, B, 2, C, 2, S, SEP, 1, M, WORK.asMatrix(2),
        2, IWORK, INFO);
    chkxer('DTRSNA', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 16;
    dtrsna('B', 'A', SEL, 2, A, 2, B, 2, C, 2, S, SEP, 2, M, WORK.asMatrix(1),
        1, IWORK, INFO);
    chkxer('DTRSNA', infoc.INFOT, NOUT, LERR, OK);
    NT += 9;
  });

  test('DTRSEN', () {
    SEL[1] = false;
    srnamc.SRNAMT = 'DTRSEN';
    infoc.INFOT = 1;
    dtrsen('X', 'N', SEL, 0, A, 1, B, 1, WR, WI, M, S.box(1), SEP.box(1), WORK,
        1, IWORK, 1, INFO);
    chkxer('DTRSEN', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 2;
    dtrsen('N', 'X', SEL, 0, A, 1, B, 1, WR, WI, M, S.box(1), SEP.box(1), WORK,
        1, IWORK, 1, INFO);
    chkxer('DTRSEN', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 4;
    dtrsen('N', 'N', SEL, -1, A, 1, B, 1, WR, WI, M, S.box(1), SEP.box(1), WORK,
        1, IWORK, 1, INFO);
    chkxer('DTRSEN', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 6;
    dtrsen('N', 'N', SEL, 2, A, 1, B, 1, WR, WI, M, S.box(1), SEP.box(1), WORK,
        2, IWORK, 1, INFO);
    chkxer('DTRSEN', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 8;
    dtrsen('N', 'V', SEL, 2, A, 2, B, 1, WR, WI, M, S.box(1), SEP.box(1), WORK,
        1, IWORK, 1, INFO);
    chkxer('DTRSEN', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 15;
    dtrsen('N', 'V', SEL, 2, A, 2, B, 2, WR, WI, M, S.box(1), SEP.box(1), WORK,
        0, IWORK, 1, INFO);
    chkxer('DTRSEN', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 15;
    dtrsen('E', 'V', SEL, 3, A, 3, B, 3, WR, WI, M, S.box(1), SEP.box(1), WORK,
        1, IWORK, 1, INFO);
    chkxer('DTRSEN', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 15;
    dtrsen('V', 'V', SEL, 3, A, 3, B, 3, WR, WI, M, S.box(1), SEP.box(1), WORK,
        3, IWORK, 2, INFO);
    chkxer('DTRSEN', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 17;
    dtrsen('E', 'V', SEL, 2, A, 2, B, 2, WR, WI, M, S.box(1), SEP.box(1), WORK,
        1, IWORK, 0, INFO);
    chkxer('DTRSEN', infoc.INFOT, NOUT, LERR, OK);
    infoc.INFOT = 17;
    dtrsen('V', 'V', SEL, 3, A, 3, B, 3, WR, WI, M, S.box(1), SEP.box(1), WORK,
        4, IWORK, 1, INFO);
    chkxer('DTRSEN', infoc.INFOT, NOUT, LERR, OK);
    NT += 10;
  });

  // Print a summary line.

  if (infoc.OK.value) {
    infoc.NOUT.println(
        ' ${PATH.a3} routines passed the tests of the error exits (${NT.i3} tests done)');
  } else {
    infoc.NOUT.println(
        ' *** ${PATH.a3} routines failed the tests of the error exits ***');
  }
}
