// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dbdsdc.dart';
import 'package:dart_lapack/src/dbdsqr.dart';
import 'package:dart_lapack/src/dbdsvdx.dart';
import 'package:dart_lapack/src/dgebd2.dart';
import 'package:dart_lapack/src/dgebrd.dart';
import 'package:dart_lapack/src/dorgbr.dart';
import 'package:dart_lapack/src/dormbr.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/lsamen.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';

import '../test_driver.dart';
import 'chkxer.dart';
import 'common.dart';

void derrbd(final String PATH, final Nout NUNIT, final TestDriver test) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 4, LW = NMAX;
  const ZERO = 0.0, ONE = 1.0;
  int NT;
  final INFO = Box(0), NS = Box(0);
  final IQ = Matrix<int>(NMAX, NMAX);
  final IW = Array<int>(NMAX);
  final A = Matrix<double>(NMAX, NMAX),
      Q = Matrix<double>(NMAX, NMAX),
      U = Matrix<double>(NMAX, NMAX),
      V = Matrix<double>(NMAX, NMAX);
  final D = Array<double>(NMAX),
      E = Array<double>(NMAX),
      S = Array<double>(NMAX),
      TP = Array<double>(NMAX),
      TQ = Array<double>(NMAX),
      W = Array<double>(LW);

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
  }
  OK.value = true;
  NT = 0;

  // Test error exits of the SVD routines.

  if (lsamen(2, C2, 'BD')) {
    test('DGEBRD', () {
      srnamc.SRNAMT = 'DGEBRD';
      infoc.INFOT = 1;
      dgebrd(-1, 0, A, 1, D, E, TQ, TP, W, 1, INFO);
      chkxer('DGEBRD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dgebrd(0, -1, A, 1, D, E, TQ, TP, W, 1, INFO);
      chkxer('DGEBRD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dgebrd(2, 1, A, 1, D, E, TQ, TP, W, 2, INFO);
      chkxer('DGEBRD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dgebrd(2, 1, A, 2, D, E, TQ, TP, W, 1, INFO);
      chkxer('DGEBRD', infoc.INFOT, NOUT, LERR, OK);
      NT += 4;
    });

    test('DGEBD2', () {
      srnamc.SRNAMT = 'DGEBD2';
      infoc.INFOT = 1;
      dgebd2(-1, 0, A, 1, D, E, TQ, TP, W, INFO);
      chkxer('DGEBD2', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dgebd2(0, -1, A, 1, D, E, TQ, TP, W, INFO);
      chkxer('DGEBD2', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dgebd2(2, 1, A, 1, D, E, TQ, TP, W, INFO);
      chkxer('DGEBD2', infoc.INFOT, NOUT, LERR, OK);
      NT += 3;
    });

    test('DORGBR', () {
      srnamc.SRNAMT = 'DORGBR';
      infoc.INFOT = 1;
      dorgbr('/', 0, 0, 0, A, 1, TQ, W, 1, INFO);
      chkxer('DORGBR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dorgbr('Q', -1, 0, 0, A, 1, TQ, W, 1, INFO);
      chkxer('DORGBR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dorgbr('Q', 0, -1, 0, A, 1, TQ, W, 1, INFO);
      chkxer('DORGBR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dorgbr('Q', 0, 1, 0, A, 1, TQ, W, 1, INFO);
      chkxer('DORGBR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dorgbr('Q', 1, 0, 1, A, 1, TQ, W, 1, INFO);
      chkxer('DORGBR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dorgbr('P', 1, 0, 0, A, 1, TQ, W, 1, INFO);
      chkxer('DORGBR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dorgbr('P', 0, 1, 1, A, 1, TQ, W, 1, INFO);
      chkxer('DORGBR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dorgbr('Q', 0, 0, -1, A, 1, TQ, W, 1, INFO);
      chkxer('DORGBR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dorgbr('Q', 2, 1, 1, A, 1, TQ, W, 1, INFO);
      chkxer('DORGBR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dorgbr('Q', 2, 2, 1, A, 2, TQ, W, 1, INFO);
      chkxer('DORGBR', infoc.INFOT, NOUT, LERR, OK);
      NT += 10;
    });

    test('DORMBR', () {
      srnamc.SRNAMT = 'DORMBR';
      infoc.INFOT = 1;
      dormbr('/', 'L', 'T', 0, 0, 0, A, 1, TQ, U, 1, W, 1, INFO);
      chkxer('DORMBR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dormbr('Q', '/', 'T', 0, 0, 0, A, 1, TQ, U, 1, W, 1, INFO);
      chkxer('DORMBR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dormbr('Q', 'L', '/', 0, 0, 0, A, 1, TQ, U, 1, W, 1, INFO);
      chkxer('DORMBR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dormbr('Q', 'L', 'T', -1, 0, 0, A, 1, TQ, U, 1, W, 1, INFO);
      chkxer('DORMBR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dormbr('Q', 'L', 'T', 0, -1, 0, A, 1, TQ, U, 1, W, 1, INFO);
      chkxer('DORMBR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dormbr('Q', 'L', 'T', 0, 0, -1, A, 1, TQ, U, 1, W, 1, INFO);
      chkxer('DORMBR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dormbr('Q', 'L', 'T', 2, 0, 0, A, 1, TQ, U, 2, W, 1, INFO);
      chkxer('DORMBR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dormbr('Q', 'R', 'T', 0, 2, 0, A, 1, TQ, U, 1, W, 1, INFO);
      chkxer('DORMBR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dormbr('P', 'L', 'T', 2, 0, 2, A, 1, TQ, U, 2, W, 1, INFO);
      chkxer('DORMBR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dormbr('P', 'R', 'T', 0, 2, 2, A, 1, TQ, U, 1, W, 1, INFO);
      chkxer('DORMBR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 11;
      dormbr('Q', 'R', 'T', 2, 0, 0, A, 1, TQ, U, 1, W, 1, INFO);
      chkxer('DORMBR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 13;
      dormbr('Q', 'L', 'T', 0, 2, 0, A, 1, TQ, U, 1, W, 1, INFO);
      chkxer('DORMBR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 13;
      dormbr('Q', 'R', 'T', 2, 0, 0, A, 1, TQ, U, 2, W, 1, INFO);
      chkxer('DORMBR', infoc.INFOT, NOUT, LERR, OK);
      NT += 13;
    });

    test('DBDSQR', () {
      srnamc.SRNAMT = 'DBDSQR';
      infoc.INFOT = 1;
      dbdsqr('/', 0, 0, 0, 0, D, E, V, 1, U, 1, A, 1, W, INFO);
      chkxer('DBDSQR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dbdsqr('U', -1, 0, 0, 0, D, E, V, 1, U, 1, A, 1, W, INFO);
      chkxer('DBDSQR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dbdsqr('U', 0, -1, 0, 0, D, E, V, 1, U, 1, A, 1, W, INFO);
      chkxer('DBDSQR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dbdsqr('U', 0, 0, -1, 0, D, E, V, 1, U, 1, A, 1, W, INFO);
      chkxer('DBDSQR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dbdsqr('U', 0, 0, 0, -1, D, E, V, 1, U, 1, A, 1, W, INFO);
      chkxer('DBDSQR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dbdsqr('U', 2, 1, 0, 0, D, E, V, 1, U, 1, A, 1, W, INFO);
      chkxer('DBDSQR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 11;
      dbdsqr('U', 0, 0, 2, 0, D, E, V, 1, U, 1, A, 1, W, INFO);
      chkxer('DBDSQR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 13;
      dbdsqr('U', 2, 0, 0, 1, D, E, V, 1, U, 1, A, 1, W, INFO);
      chkxer('DBDSQR', infoc.INFOT, NOUT, LERR, OK);
      NT += 8;
    });

    test('DBDSDC', () {
      srnamc.SRNAMT = 'DBDSDC';
      infoc.INFOT = 1;
      dbdsdc('/', 'N', 0, D, E, U, 1, V, 1, Q.asArray(), IQ.asArray(), W, IW,
          INFO);
      chkxer('DBDSDC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dbdsdc('U', '/', 0, D, E, U, 1, V, 1, Q.asArray(), IQ.asArray(), W, IW,
          INFO);
      chkxer('DBDSDC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dbdsdc('U', 'N', -1, D, E, U, 1, V, 1, Q.asArray(), IQ.asArray(), W, IW,
          INFO);
      chkxer('DBDSDC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dbdsdc('U', 'I', 2, D, E, U, 1, V, 1, Q.asArray(), IQ.asArray(), W, IW,
          INFO);
      chkxer('DBDSDC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dbdsdc('U', 'I', 2, D, E, U, 2, V, 1, Q.asArray(), IQ.asArray(), W, IW,
          INFO);
      chkxer('DBDSDC', infoc.INFOT, NOUT, LERR, OK);
      NT += 5;
    });

    test('DBDSVDX', () {
      srnamc.SRNAMT = 'DBDSVDX';
      infoc.INFOT = 1;
      dbdsvdx(
          'X', 'N', 'A', 1, D, E, ZERO, ONE, 0, 0, NS, S, Q, 1, W, IW, INFO);
      chkxer('DBDSVDX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dbdsvdx(
          'U', 'X', 'A', 1, D, E, ZERO, ONE, 0, 0, NS, S, Q, 1, W, IW, INFO);
      chkxer('DBDSVDX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dbdsvdx(
          'U', 'V', 'X', 1, D, E, ZERO, ONE, 0, 0, NS, S, Q, 1, W, IW, INFO);
      chkxer('DBDSVDX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dbdsvdx(
          'U', 'V', 'A', -1, D, E, ZERO, ONE, 0, 0, NS, S, Q, 1, W, IW, INFO);
      chkxer('DBDSVDX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dbdsvdx(
          'U', 'V', 'V', 2, D, E, -ONE, ZERO, 0, 0, NS, S, Q, 1, W, IW, INFO);
      chkxer('DBDSVDX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dbdsvdx(
          'U', 'V', 'V', 2, D, E, ONE, ZERO, 0, 0, NS, S, Q, 1, W, IW, INFO);
      chkxer('DBDSVDX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dbdsvdx(
          'L', 'V', 'I', 2, D, E, ZERO, ZERO, 0, 2, NS, S, Q, 1, W, IW, INFO);
      chkxer('DBDSVDX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dbdsvdx(
          'L', 'V', 'I', 4, D, E, ZERO, ZERO, 5, 2, NS, S, Q, 1, W, IW, INFO);
      chkxer('DBDSVDX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dbdsvdx(
          'L', 'V', 'I', 4, D, E, ZERO, ZERO, 3, 2, NS, S, Q, 1, W, IW, INFO);
      chkxer('DBDSVDX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dbdsvdx(
          'L', 'V', 'I', 4, D, E, ZERO, ZERO, 3, 5, NS, S, Q, 1, W, IW, INFO);
      chkxer('DBDSVDX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 14;
      dbdsvdx(
          'L', 'V', 'A', 4, D, E, ZERO, ZERO, 0, 0, NS, S, Q, 0, W, IW, INFO);
      chkxer('DBDSVDX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 14;
      dbdsvdx(
          'L', 'V', 'A', 4, D, E, ZERO, ZERO, 0, 0, NS, S, Q, 2, W, IW, INFO);
      chkxer('DBDSVDX', infoc.INFOT, NOUT, LERR, OK);
      NT += 12;
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
