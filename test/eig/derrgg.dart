// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgges.dart';
import 'package:lapack/src/dgges3.dart';
import 'package:lapack/src/dggesx.dart';
import 'package:lapack/src/dggev.dart';
import 'package:lapack/src/dggev3.dart';
import 'package:lapack/src/dggevx.dart';
import 'package:lapack/src/dggglm.dart';
import 'package:lapack/src/dgghd3.dart';
import 'package:lapack/src/dgghrd.dart';
import 'package:lapack/src/dgglse.dart';
import 'package:lapack/src/dggqrf.dart';
import 'package:lapack/src/dggrqf.dart';
import 'package:lapack/src/dggsvd3.dart';
import 'package:lapack/src/dggsvp3.dart';
import 'package:lapack/src/dhgeqz.dart';
import 'package:lapack/src/dorcsd.dart';
import 'package:lapack/src/dtgevc.dart';
import 'package:lapack/src/dtgexc.dart';
import 'package:lapack/src/dtgsen.dart';
import 'package:lapack/src/dtgsja.dart';
import 'package:lapack/src/dtgsna.dart';
import 'package:lapack/src/dtgsyl.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../test_driver.dart';
import 'chkxer.dart';
import 'common.dart';
import 'dlctes.dart';
import 'dlctsx.dart';
import 'xlaenv.dart';

void derrgg(final String PATH, final Nout NUNIT, final TestDriver test) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 3, LW = 6 * NMAX;
  const ONE = 1.0, ZERO = 0.0;
  int NT;
  final BW = Array<bool>(NMAX), SEL = Array<bool>(NMAX);
  final IW = Array<int>(NMAX), IDUM = Array<int>(NMAX);
  final A = Matrix<double>(NMAX, NMAX),
      B = Matrix<double>(NMAX, NMAX),
      Q = Matrix<double>(NMAX, NMAX),
      U = Matrix<double>(NMAX, NMAX),
      V = Matrix<double>(NMAX, NMAX),
      Z = Matrix<double>(NMAX, NMAX);
  final LS = Array<double>(NMAX),
      R1 = Array<double>(NMAX),
      R2 = Array<double>(NMAX),
      R3 = Array<double>(NMAX),
      RCE = Array<double>(2),
      RCV = Array<double>(2),
      RS = Array<double>(NMAX),
      TAU = Array<double>(NMAX),
      W = Array<double>(LW);
  final INFO = Box(0),
      M = Box(0),
      DUMMYK = Box(0),
      DUMMYL = Box(0),
      NCYCLE = Box(0),
      SDIM = Box(0),
      ILO = Box(0),
      IHI = Box(0),
      IFST = Box(0),
      ILST = Box(0);
  final DIF = Box(0.0),
      SCALE = Box(0.0),
      ANRM = Box(0.0),
      BNRM = Box(0.0),
      TOLA = Box(0.0),
      TOLB = Box(0.0);

  final NOUT = infoc.NOUT = NUNIT;
  final OK = infoc.OK;
  final LERR = infoc.LERR;
  NOUT.println();
  final C2 = PATH.substring(1, 3);

  // Set the variables to innocuous values.

  for (var J = 1; J <= NMAX; J++) {
    SEL[J] = true;
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = ZERO;
      B[I][J] = ZERO;
    }
  }
  for (var I = 1; I <= NMAX; I++) {
    A[I][I] = ONE;
    B[I][I] = ONE;
  }
  infoc.OK.value = true;
  TOLA.value = 1.0;
  TOLB.value = 1.0;
  IFST.value = 1;
  ILST.value = 1;
  NT = 0;
  const LWORK = 1;

  // Call XLAENV to set the parameters used in CLAQZ0

  xlaenv(12, 10);
  xlaenv(13, 12);
  xlaenv(14, 13);
  xlaenv(15, 2);
  xlaenv(17, 10);

  // Test error exits for the GG path.

  if (lsamen(2, C2, 'GG')) {
    test('DGGHRD', () {
      srnamc.SRNAMT = 'DGGHRD';
      infoc.INFOT = 1;
      dgghrd('/', 'N', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, INFO);
      chkxer('DGGHRD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dgghrd('N', '/', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, INFO);
      chkxer('DGGHRD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dgghrd('N', 'N', -1, 0, 0, A, 1, B, 1, Q, 1, Z, 1, INFO);
      chkxer('DGGHRD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dgghrd('N', 'N', 0, 0, 0, A, 1, B, 1, Q, 1, Z, 1, INFO);
      chkxer('DGGHRD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dgghrd('N', 'N', 0, 1, 1, A, 1, B, 1, Q, 1, Z, 1, INFO);
      chkxer('DGGHRD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dgghrd('N', 'N', 2, 1, 1, A, 1, B, 2, Q, 1, Z, 1, INFO);
      chkxer('DGGHRD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dgghrd('N', 'N', 2, 1, 1, A, 2, B, 1, Q, 1, Z, 1, INFO);
      chkxer('DGGHRD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 11;
      dgghrd('V', 'N', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, INFO);
      chkxer('DGGHRD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 13;
      dgghrd('N', 'V', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, INFO);
      chkxer('DGGHRD', infoc.INFOT, NOUT, LERR, OK);
      NT += 9;
    });

    test('DGGHD3', () {
      srnamc.SRNAMT = 'DGGHD3';
      infoc.INFOT = 1;
      dgghd3('/', 'N', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO);
      chkxer('DGGHD3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dgghd3('N', '/', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO);
      chkxer('DGGHD3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dgghd3('N', 'N', -1, 0, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO);
      chkxer('DGGHD3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dgghd3('N', 'N', 0, 0, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO);
      chkxer('DGGHD3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dgghd3('N', 'N', 0, 1, 1, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO);
      chkxer('DGGHD3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dgghd3('N', 'N', 2, 1, 1, A, 1, B, 2, Q, 1, Z, 1, W, LW, INFO);
      chkxer('DGGHD3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dgghd3('N', 'N', 2, 1, 1, A, 2, B, 1, Q, 1, Z, 1, W, LW, INFO);
      chkxer('DGGHD3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 11;
      dgghd3('V', 'N', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, W, LW, INFO);
      chkxer('DGGHD3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 13;
      dgghd3('N', 'V', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, W, LW, INFO);
      chkxer('DGGHD3', infoc.INFOT, NOUT, LERR, OK);
      NT += 9;
    });

    test('DHGEQZ', () {
      srnamc.SRNAMT = 'DHGEQZ';
      infoc.INFOT = 1;
      dhgeqz('/', 'N', 'N', 0, 1, 0, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, W, LW,
          INFO);
      chkxer('DHGEQZ', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dhgeqz('E', '/', 'N', 0, 1, 0, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, W, LW,
          INFO);
      chkxer('DHGEQZ', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dhgeqz('E', 'N', '/', 0, 1, 0, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, W, LW,
          INFO);
      chkxer('DHGEQZ', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dhgeqz('E', 'N', 'N', -1, 0, 0, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, W, LW,
          INFO);
      chkxer('DHGEQZ', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dhgeqz('E', 'N', 'N', 0, 0, 0, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, W, LW,
          INFO);
      chkxer('DHGEQZ', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dhgeqz('E', 'N', 'N', 0, 1, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, W, LW,
          INFO);
      chkxer('DHGEQZ', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dhgeqz('E', 'N', 'N', 2, 1, 1, A, 1, B, 2, R1, R2, R3, Q, 1, Z, 1, W, LW,
          INFO);
      chkxer('DHGEQZ', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dhgeqz('E', 'N', 'N', 2, 1, 1, A, 2, B, 1, R1, R2, R3, Q, 1, Z, 1, W, LW,
          INFO);
      chkxer('DHGEQZ', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 15;
      dhgeqz('E', 'V', 'N', 2, 1, 1, A, 2, B, 2, R1, R2, R3, Q, 1, Z, 1, W, LW,
          INFO);
      chkxer('DHGEQZ', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 17;
      dhgeqz('E', 'N', 'V', 2, 1, 1, A, 2, B, 2, R1, R2, R3, Q, 1, Z, 1, W, LW,
          INFO);
      chkxer('DHGEQZ', infoc.INFOT, NOUT, LERR, OK);
      NT += 10;
    });

    test('DTGEVC', () {
      srnamc.SRNAMT = 'DTGEVC';
      infoc.INFOT = 1;
      dtgevc('/', 'A', SEL, 0, A, 1, B, 1, Q, 1, Z, 1, 0, M, W, INFO);
      chkxer('DTGEVC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dtgevc('R', '/', SEL, 0, A, 1, B, 1, Q, 1, Z, 1, 0, M, W, INFO);
      chkxer('DTGEVC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dtgevc('R', 'A', SEL, -1, A, 1, B, 1, Q, 1, Z, 1, 0, M, W, INFO);
      chkxer('DTGEVC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dtgevc('R', 'A', SEL, 2, A, 1, B, 2, Q, 1, Z, 2, 0, M, W, INFO);
      chkxer('DTGEVC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dtgevc('R', 'A', SEL, 2, A, 2, B, 1, Q, 1, Z, 2, 0, M, W, INFO);
      chkxer('DTGEVC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dtgevc('L', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 1, 0, M, W, INFO);
      chkxer('DTGEVC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 12;
      dtgevc('R', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 1, 0, M, W, INFO);
      chkxer('DTGEVC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 13;
      dtgevc('R', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 2, 1, M, W, INFO);
      chkxer('DTGEVC', infoc.INFOT, NOUT, LERR, OK);
      NT += 8;
    });

    // Test error exits for the GSV path.
  } else if (lsamen(3, PATH, 'GSV')) {
    test('DGGSVD3', () {
      srnamc.SRNAMT = 'DGGSVD3';
      infoc.INFOT = 1;
      dggsvd3('/', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1,
          V, 1, Q, 1, W, LWORK, IDUM, INFO);
      chkxer('DGGSVD3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dggsvd3('N', '/', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1,
          V, 1, Q, 1, W, LWORK, IDUM, INFO);
      chkxer('DGGSVD3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dggsvd3('N', 'N', '/', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1,
          V, 1, Q, 1, W, LWORK, IDUM, INFO);
      chkxer('DGGSVD3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dggsvd3('N', 'N', 'N', -1, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1,
          V, 1, Q, 1, W, LWORK, IDUM, INFO);
      chkxer('DGGSVD3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dggsvd3('N', 'N', 'N', 0, -1, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1,
          V, 1, Q, 1, W, LWORK, IDUM, INFO);
      chkxer('DGGSVD3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dggsvd3('N', 'N', 'N', 0, 0, -1, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1,
          V, 1, Q, 1, W, LWORK, IDUM, INFO);
      chkxer('DGGSVD3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dggsvd3('N', 'N', 'N', 2, 1, 1, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1,
          V, 1, Q, 1, W, LWORK, IDUM, INFO);
      chkxer('DGGSVD3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 12;
      dggsvd3('N', 'N', 'N', 1, 1, 2, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1,
          V, 1, Q, 1, W, LWORK, IDUM, INFO);
      chkxer('DGGSVD3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 16;
      dggsvd3('U', 'N', 'N', 2, 2, 2, DUMMYK, DUMMYL, A, 2, B, 2, R1, R2, U, 1,
          V, 1, Q, 1, W, LWORK, IDUM, INFO);
      chkxer('DGGSVD3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 18;
      dggsvd3('N', 'V', 'N', 1, 1, 2, DUMMYK, DUMMYL, A, 1, B, 2, R1, R2, U, 1,
          V, 1, Q, 1, W, LWORK, IDUM, INFO);
      chkxer('DGGSVD3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 20;
      dggsvd3('N', 'N', 'Q', 1, 2, 1, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1,
          V, 1, Q, 1, W, LWORK, IDUM, INFO);
      chkxer('DGGSVD3', infoc.INFOT, NOUT, LERR, OK);
      NT += 11;
    });

    test('DGGSVP3', () {
      srnamc.SRNAMT = 'DGGSVP3';
      infoc.INFOT = 1;
      dggsvp3('/', 'N', 'N', 0, 0, 0, A, 1, B, 1, TOLA.value, TOLB.value,
          DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO);
      chkxer('DGGSVP3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dggsvp3('N', '/', 'N', 0, 0, 0, A, 1, B, 1, TOLA.value, TOLB.value,
          DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO);
      chkxer('DGGSVP3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dggsvp3('N', 'N', '/', 0, 0, 0, A, 1, B, 1, TOLA.value, TOLB.value,
          DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO);
      chkxer('DGGSVP3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dggsvp3('N', 'N', 'N', -1, 0, 0, A, 1, B, 1, TOLA.value, TOLB.value,
          DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO);
      chkxer('DGGSVP3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dggsvp3('N', 'N', 'N', 0, -1, 0, A, 1, B, 1, TOLA.value, TOLB.value,
          DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO);
      chkxer('DGGSVP3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dggsvp3('N', 'N', 'N', 0, 0, -1, A, 1, B, 1, TOLA.value, TOLB.value,
          DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO);
      chkxer('DGGSVP3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dggsvp3('N', 'N', 'N', 2, 1, 1, A, 1, B, 1, TOLA.value, TOLB.value,
          DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO);
      chkxer('DGGSVP3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dggsvp3('N', 'N', 'N', 1, 2, 1, A, 1, B, 1, TOLA.value, TOLB.value,
          DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO);
      chkxer('DGGSVP3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 16;
      dggsvp3('U', 'N', 'N', 2, 2, 2, A, 2, B, 2, TOLA.value, TOLB.value,
          DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO);
      chkxer('DGGSVP3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 18;
      dggsvp3('N', 'V', 'N', 1, 2, 1, A, 1, B, 2, TOLA.value, TOLB.value,
          DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO);
      chkxer('DGGSVP3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 20;
      dggsvp3('N', 'N', 'Q', 1, 1, 2, A, 1, B, 1, TOLA.value, TOLB.value,
          DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO);
      chkxer('DGGSVP3', infoc.INFOT, NOUT, LERR, OK);
      NT += 11;
    });

    test('DTGSJA', () {
      srnamc.SRNAMT = 'DTGSJA';
      infoc.INFOT = 1;
      dtgsja('/', 'N', 'N', 0, 0, 0, DUMMYK.value, DUMMYL.value, A, 1, B, 1,
          TOLA.value, TOLB.value, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO);
      chkxer('DTGSJA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dtgsja('N', '/', 'N', 0, 0, 0, DUMMYK.value, DUMMYL.value, A, 1, B, 1,
          TOLA.value, TOLB.value, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO);
      chkxer('DTGSJA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dtgsja('N', 'N', '/', 0, 0, 0, DUMMYK.value, DUMMYL.value, A, 1, B, 1,
          TOLA.value, TOLB.value, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO);
      chkxer('DTGSJA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dtgsja('N', 'N', 'N', -1, 0, 0, DUMMYK.value, DUMMYL.value, A, 1, B, 1,
          TOLA.value, TOLB.value, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO);
      chkxer('DTGSJA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dtgsja('N', 'N', 'N', 0, -1, 0, DUMMYK.value, DUMMYL.value, A, 1, B, 1,
          TOLA.value, TOLB.value, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO);
      chkxer('DTGSJA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dtgsja('N', 'N', 'N', 0, 0, -1, DUMMYK.value, DUMMYL.value, A, 1, B, 1,
          TOLA.value, TOLB.value, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO);
      chkxer('DTGSJA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dtgsja('N', 'N', 'N', 0, 0, 0, DUMMYK.value, DUMMYL.value, A, 0, B, 1,
          TOLA.value, TOLB.value, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO);
      chkxer('DTGSJA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 12;
      dtgsja('N', 'N', 'N', 0, 0, 0, DUMMYK.value, DUMMYL.value, A, 1, B, 0,
          TOLA.value, TOLB.value, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO);
      chkxer('DTGSJA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 18;
      dtgsja('U', 'N', 'N', 0, 0, 0, DUMMYK.value, DUMMYL.value, A, 1, B, 1,
          TOLA.value, TOLB.value, R1, R2, U, 0, V, 1, Q, 1, W, NCYCLE, INFO);
      chkxer('DTGSJA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 20;
      dtgsja('N', 'V', 'N', 0, 0, 0, DUMMYK.value, DUMMYL.value, A, 1, B, 1,
          TOLA.value, TOLB.value, R1, R2, U, 1, V, 0, Q, 1, W, NCYCLE, INFO);
      chkxer('DTGSJA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 22;
      dtgsja('N', 'N', 'Q', 0, 0, 0, DUMMYK.value, DUMMYL.value, A, 1, B, 1,
          TOLA.value, TOLB.value, R1, R2, U, 1, V, 1, Q, 0, W, NCYCLE, INFO);
      chkxer('DTGSJA', infoc.INFOT, NOUT, LERR, OK);
      NT += 11;
    });

    // Test error exits for the GLM path.
  } else if (lsamen(3, PATH, 'GLM')) {
    test('DGGGLM', () {
      srnamc.SRNAMT = 'DGGGLM';
      infoc.INFOT = 1;
      dggglm(-1, 0, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO);
      chkxer('DGGGLM', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dggglm(0, -1, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO);
      chkxer('DGGGLM', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dggglm(0, 1, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO);
      chkxer('DGGGLM', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dggglm(0, 0, -1, A, 1, B, 1, R1, R2, R3, W, LW, INFO);
      chkxer('DGGGLM', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dggglm(1, 0, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO);
      chkxer('DGGGLM', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dggglm(0, 0, 0, A, 0, B, 1, R1, R2, R3, W, LW, INFO);
      chkxer('DGGGLM', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dggglm(0, 0, 0, A, 1, B, 0, R1, R2, R3, W, LW, INFO);
      chkxer('DGGGLM', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 12;
      dggglm(1, 1, 1, A, 1, B, 1, R1, R2, R3, W, 1, INFO);
      chkxer('DGGGLM', infoc.INFOT, NOUT, LERR, OK);
      NT += 8;
    });

    // Test error exits for the LSE path.
  } else if (lsamen(3, PATH, 'LSE')) {
    test('DGGLSE', () {
      srnamc.SRNAMT = 'DGGLSE';
      infoc.INFOT = 1;
      dgglse(-1, 0, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO);
      chkxer('DGGLSE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dgglse(0, -1, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO);
      chkxer('DGGLSE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dgglse(0, 0, -1, A, 1, B, 1, R1, R2, R3, W, LW, INFO);
      chkxer('DGGLSE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dgglse(0, 0, 1, A, 1, B, 1, R1, R2, R3, W, LW, INFO);
      chkxer('DGGLSE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dgglse(0, 1, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO);
      chkxer('DGGLSE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dgglse(0, 0, 0, A, 0, B, 1, R1, R2, R3, W, LW, INFO);
      chkxer('DGGLSE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dgglse(0, 0, 0, A, 1, B, 0, R1, R2, R3, W, LW, INFO);
      chkxer('DGGLSE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 12;
      dgglse(1, 1, 1, A, 1, B, 1, R1, R2, R3, W, 1, INFO);
      chkxer('DGGLSE', infoc.INFOT, NOUT, LERR, OK);
      NT += 8;
    });

    // Test error exits for the CSD path.
  } else if (lsamen(3, PATH, 'CSD')) {
    test('DORCSD', () {
      srnamc.SRNAMT = 'DORCSD';
      infoc.INFOT = 7;
      dorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', -1, 0, 0, A, 1, A, 1, A, 1, A, 1,
          A.asArray(), A, 1, A, 1, A, 1, A, 1, W, LW, IW, INFO);
      chkxer('DORCSD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, -1, 0, A, 1, A, 1, A, 1, A, 1,
          A.asArray(), A, 1, A, 1, A, 1, A, 1, W, LW, IW, INFO);
      chkxer('DORCSD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, -1, A, 1, A, 1, A, 1, A, 1,
          A.asArray(), A, 1, A, 1, A, 1, A, 1, W, LW, IW, INFO);
      chkxer('DORCSD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 11;
      dorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, -1, A, 1, A, 1, A, 1,
          A.asArray(), A, 1, A, 1, A, 1, A, 1, W, LW, IW, INFO);
      chkxer('DORCSD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 20;
      dorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1,
          A.asArray(), A, -1, A, 1, A, 1, A, 1, W, LW, IW, INFO);
      chkxer('DORCSD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 22;
      dorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1,
          A.asArray(), A, 1, A, -1, A, 1, A, 1, W, LW, IW, INFO);
      chkxer('DORCSD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 24;
      dorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1,
          A.asArray(), A, 1, A, 1, A, -1, A, 1, W, LW, IW, INFO);
      chkxer('DORCSD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 26;
      dorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1,
          A.asArray(), A, 1, A, 1, A, 1, A, -1, W, LW, IW, INFO);
      chkxer('DORCSD', infoc.INFOT, NOUT, LERR, OK);
      NT += 8;
    });

    // Test error exits for the GQR path.
  } else if (lsamen(3, PATH, 'GQR')) {
    test('DGGQRF', () {
      srnamc.SRNAMT = 'DGGQRF';
      infoc.INFOT = 1;
      dggqrf(-1, 0, 0, A, 1, R1, B, 1, R2, W, LW, INFO);
      chkxer('DGGQRF', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dggqrf(0, -1, 0, A, 1, R1, B, 1, R2, W, LW, INFO);
      chkxer('DGGQRF', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dggqrf(0, 0, -1, A, 1, R1, B, 1, R2, W, LW, INFO);
      chkxer('DGGQRF', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dggqrf(0, 0, 0, A, 0, R1, B, 1, R2, W, LW, INFO);
      chkxer('DGGQRF', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dggqrf(0, 0, 0, A, 1, R1, B, 0, R2, W, LW, INFO);
      chkxer('DGGQRF', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 11;
      dggqrf(1, 1, 2, A, 1, R1, B, 1, R2, W, 1, INFO);
      chkxer('DGGQRF', infoc.INFOT, NOUT, LERR, OK);
      NT += 6;
    });

    test('DGGRQF', () {
      srnamc.SRNAMT = 'DGGRQF';
      infoc.INFOT = 1;
      dggrqf(-1, 0, 0, A, 1, R1, B, 1, R2, W, LW, INFO);
      chkxer('DGGRQF', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dggrqf(0, -1, 0, A, 1, R1, B, 1, R2, W, LW, INFO);
      chkxer('DGGRQF', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dggrqf(0, 0, -1, A, 1, R1, B, 1, R2, W, LW, INFO);
      chkxer('DGGRQF', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dggrqf(0, 0, 0, A, 0, R1, B, 1, R2, W, LW, INFO);
      chkxer('DGGRQF', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dggrqf(0, 0, 0, A, 1, R1, B, 0, R2, W, LW, INFO);
      chkxer('DGGRQF', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 11;
      dggrqf(1, 1, 2, A, 1, R1, B, 1, R2, W, 1, INFO);
      chkxer('DGGRQF', infoc.INFOT, NOUT, LERR, OK);
      NT += 6;
    });

    // Test error exits for the DGS, DGV, DGX, and DXV paths.
  } else if (lsamen(3, PATH, 'DGS') ||
      lsamen(3, PATH, 'DGV') ||
      lsamen(3, PATH, 'DGX') ||
      lsamen(3, PATH, 'DXV')) {
    test('DGGES', () {
      srnamc.SRNAMT = 'DGGES';
      infoc.INFOT = 1;
      dgges('/', 'N', 'S', dlctes, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1,
          W, 1, BW, INFO);
      chkxer('DGGES', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dgges('N', '/', 'S', dlctes, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1,
          W, 1, BW, INFO);
      chkxer('DGGES', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dgges('N', 'V', '/', dlctes, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1,
          W, 1, BW, INFO);
      chkxer('DGGES', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dgges('N', 'V', 'S', dlctes, -1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1,
          W, 1, BW, INFO);
      chkxer('DGGES', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dgges('N', 'V', 'S', dlctes, 1, A, 0, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1,
          W, 1, BW, INFO);
      chkxer('DGGES', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dgges('N', 'V', 'S', dlctes, 1, A, 1, B, 0, SDIM, R1, R2, R3, Q, 1, U, 1,
          W, 1, BW, INFO);
      chkxer('DGGES', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 15;
      dgges('N', 'V', 'S', dlctes, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 0, U, 1,
          W, 1, BW, INFO);
      chkxer('DGGES', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 15;
      dgges('V', 'V', 'S', dlctes, 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 1, U, 2,
          W, 1, BW, INFO);
      chkxer('DGGES', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 17;
      dgges('N', 'V', 'S', dlctes, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 0,
          W, 1, BW, INFO);
      chkxer('DGGES', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 17;
      dgges('V', 'V', 'S', dlctes, 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 2, U, 1,
          W, 1, BW, INFO);
      chkxer('DGGES', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 19;
      dgges('V', 'V', 'S', dlctes, 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 2, U, 2,
          W, 1, BW, INFO);
      chkxer('DGGES', infoc.INFOT, NOUT, LERR, OK);
      NT += 11;
    });

    test('DGGES3', () {
      srnamc.SRNAMT = 'DGGES3';
      infoc.INFOT = 1;
      dgges3('/', 'N', 'S', dlctes, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1,
          W, 1, BW, INFO);
      chkxer('DGGES3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dgges3('N', '/', 'S', dlctes, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1,
          W, 1, BW, INFO);
      chkxer('DGGES3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dgges3('N', 'V', '/', dlctes, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1,
          W, 1, BW, INFO);
      chkxer('DGGES3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dgges3('N', 'V', 'S', dlctes, -1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U,
          1, W, 1, BW, INFO);
      chkxer('DGGES3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dgges3('N', 'V', 'S', dlctes, 1, A, 0, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1,
          W, 1, BW, INFO);
      chkxer('DGGES3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dgges3('N', 'V', 'S', dlctes, 1, A, 1, B, 0, SDIM, R1, R2, R3, Q, 1, U, 1,
          W, 1, BW, INFO);
      chkxer('DGGES3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 15;
      dgges3('N', 'V', 'S', dlctes, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 0, U, 1,
          W, 1, BW, INFO);
      chkxer('DGGES3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 15;
      dgges3('V', 'V', 'S', dlctes, 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 1, U, 2,
          W, 1, BW, INFO);
      chkxer('DGGES3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 17;
      dgges3('N', 'V', 'S', dlctes, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 0,
          W, 1, BW, INFO);
      chkxer('DGGES3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 17;
      dgges3('V', 'V', 'S', dlctes, 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 2, U, 1,
          W, 1, BW, INFO);
      chkxer('DGGES3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 19;
      dgges3('V', 'V', 'S', dlctes, 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 2, U, 2,
          W, 1, BW, INFO);
      chkxer('DGGES3', infoc.INFOT, NOUT, LERR, OK);
      NT += 11;
    });

    test('DGGESX', () {
      srnamc.SRNAMT = 'DGGESX';
      infoc.INFOT = 1;
      dggesx('/', 'N', 'S', dlctsx, 'N', 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1,
          U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO);
      chkxer('DGGESX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dggesx('N', '/', 'S', dlctsx, 'N', 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1,
          U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO);
      chkxer('DGGESX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dggesx('V', 'V', '/', dlctsx, 'N', 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1,
          U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO);
      chkxer('DGGESX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dggesx('V', 'V', 'S', dlctsx, '/', 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1,
          U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO);
      chkxer('DGGESX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dggesx('V', 'V', 'S', dlctsx, 'B', -1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1,
          U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO);
      chkxer('DGGESX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dggesx('V', 'V', 'S', dlctsx, 'B', 1, A, 0, B, 1, SDIM, R1, R2, R3, Q, 1,
          U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO);
      chkxer('DGGESX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dggesx('V', 'V', 'S', dlctsx, 'B', 1, A, 1, B, 0, SDIM, R1, R2, R3, Q, 1,
          U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO);
      chkxer('DGGESX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 16;
      dggesx('V', 'V', 'S', dlctsx, 'B', 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 0,
          U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO);
      chkxer('DGGESX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 16;
      dggesx('V', 'V', 'S', dlctsx, 'B', 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 1,
          U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO);
      chkxer('DGGESX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 18;
      dggesx('V', 'V', 'S', dlctsx, 'B', 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1,
          U, 0, RCE, RCV, W, 1, IW, 1, BW, INFO);
      chkxer('DGGESX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 18;
      dggesx('V', 'V', 'S', dlctsx, 'B', 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 2,
          U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO);
      chkxer('DGGESX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 22;
      dggesx('V', 'V', 'S', dlctsx, 'B', 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 2,
          U, 2, RCE, RCV, W, 1, IW, 1, BW, INFO);
      chkxer('DGGESX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 24;
      dggesx('V', 'V', 'S', dlctsx, 'V', 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1,
          U, 1, RCE, RCV, W, 32, IW, 0, BW, INFO);
      chkxer('DGGESX', infoc.INFOT, NOUT, LERR, OK);
      NT += 13;
    });

    test('DGGEV', () {
      srnamc.SRNAMT = 'DGGEV';
      infoc.INFOT = 1;
      dggev('/', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO);
      chkxer('DGGEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dggev('N', '/', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO);
      chkxer('DGGEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dggev('V', 'V', -1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO);
      chkxer('DGGEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dggev('V', 'V', 1, A, 0, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO);
      chkxer('DGGEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dggev('V', 'V', 1, A, 1, B, 0, R1, R2, R3, Q, 1, U, 1, W, 1, INFO);
      chkxer('DGGEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 12;
      dggev('N', 'V', 1, A, 1, B, 1, R1, R2, R3, Q, 0, U, 1, W, 1, INFO);
      chkxer('DGGEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 12;
      dggev('V', 'V', 2, A, 2, B, 2, R1, R2, R3, Q, 1, U, 2, W, 1, INFO);
      chkxer('DGGEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 14;
      dggev('V', 'N', 2, A, 2, B, 2, R1, R2, R3, Q, 2, U, 0, W, 1, INFO);
      chkxer('DGGEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 14;
      dggev('V', 'V', 2, A, 2, B, 2, R1, R2, R3, Q, 2, U, 1, W, 1, INFO);
      chkxer('DGGEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 16;
      dggev('V', 'V', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO);
      chkxer('DGGEV', infoc.INFOT, NOUT, LERR, OK);
      NT += 10;
    });

    test('DGGEV3', () {
      xlaenv(12, 20);
      xlaenv(13, 4);
      xlaenv(14, 13);
      xlaenv(15, 2);
      xlaenv(17, 10);
      srnamc.SRNAMT = 'DGGEV3';
      infoc.INFOT = 1;
      dggev3('/', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO);
      chkxer('DGGEV3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dggev3('N', '/', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO);
      chkxer('DGGEV3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dggev3('V', 'V', -1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO);
      chkxer('DGGEV3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dggev3('V', 'V', 1, A, 0, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO);
      chkxer('DGGEV3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dggev3('V', 'V', 1, A, 1, B, 0, R1, R2, R3, Q, 1, U, 1, W, 1, INFO);
      chkxer('DGGEV3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 12;
      dggev3('N', 'V', 1, A, 1, B, 1, R1, R2, R3, Q, 0, U, 1, W, 1, INFO);
      chkxer('DGGEV3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 12;
      dggev3('V', 'V', 2, A, 2, B, 2, R1, R2, R3, Q, 1, U, 2, W, 1, INFO);
      chkxer('DGGEV3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 14;
      dggev3('V', 'N', 2, A, 2, B, 2, R1, R2, R3, Q, 2, U, 0, W, 1, INFO);
      chkxer('DGGEV3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 14;
      dggev3('V', 'V', 2, A, 2, B, 2, R1, R2, R3, Q, 2, U, 1, W, 1, INFO);
      chkxer('DGGEV3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 16;
      dggev3('V', 'V', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO);
      chkxer('DGGEV3', infoc.INFOT, NOUT, LERR, OK);
      NT += 10;
    });

    test('DGGEVX', () {
      srnamc.SRNAMT = 'DGGEVX';
      infoc.INFOT = 1;
      dggevx('/', 'N', 'N', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, ILO,
          IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO);
      chkxer('DGGEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dggevx('N', '/', 'N', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, ILO,
          IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO);
      chkxer('DGGEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dggevx('N', 'N', '/', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, ILO,
          IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO);
      chkxer('DGGEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dggevx('N', 'N', 'N', '/', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, ILO,
          IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO);
      chkxer('DGGEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dggevx('N', 'N', 'N', 'N', -1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, ILO,
          IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO);
      chkxer('DGGEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dggevx('N', 'N', 'N', 'N', 1, A, 0, B, 1, R1, R2, R3, Q, 1, U, 1, ILO,
          IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO);
      chkxer('DGGEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dggevx('N', 'N', 'N', 'N', 1, A, 1, B, 0, R1, R2, R3, Q, 1, U, 1, ILO,
          IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO);
      chkxer('DGGEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 14;
      dggevx('N', 'N', 'N', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 0, U, 1, ILO,
          IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO);
      chkxer('DGGEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 14;
      dggevx('N', 'V', 'N', 'N', 2, A, 2, B, 2, R1, R2, R3, Q, 1, U, 2, ILO,
          IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO);
      chkxer('DGGEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 16;
      dggevx('N', 'N', 'N', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 0, ILO,
          IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO);
      chkxer('DGGEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 16;
      dggevx('N', 'N', 'V', 'N', 2, A, 2, B, 2, R1, R2, R3, Q, 2, U, 1, ILO,
          IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO);
      chkxer('DGGEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 26;
      dggevx('N', 'N', 'V', 'N', 2, A, 2, B, 2, R1, R2, R3, Q, 2, U, 2, ILO,
          IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO);
      chkxer('DGGEVX', infoc.INFOT, NOUT, LERR, OK);
      NT += 12;
    });

    test('DTGEXC', () {
      srnamc.SRNAMT = 'DTGEXC';
      infoc.INFOT = 3;
      dtgexc(true, true, -1, A, 1, B, 1, Q, 1, Z, 1, IFST, ILST, W, 1, INFO);
      chkxer('DTGEXC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dtgexc(true, true, 1, A, 0, B, 1, Q, 1, Z, 1, IFST, ILST, W, 1, INFO);
      chkxer('DTGEXC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dtgexc(true, true, 1, A, 1, B, 0, Q, 1, Z, 1, IFST, ILST, W, 1, INFO);
      chkxer('DTGEXC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dtgexc(false, true, 1, A, 1, B, 1, Q, 0, Z, 1, IFST, ILST, W, 1, INFO);
      chkxer('DTGEXC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dtgexc(true, true, 1, A, 1, B, 1, Q, 0, Z, 1, IFST, ILST, W, 1, INFO);
      chkxer('DTGEXC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 11;
      dtgexc(true, false, 1, A, 1, B, 1, Q, 1, Z, 0, IFST, ILST, W, 1, INFO);
      chkxer('DTGEXC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 11;
      dtgexc(true, true, 1, A, 1, B, 1, Q, 1, Z, 0, IFST, ILST, W, 1, INFO);
      chkxer('DTGEXC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 15;
      dtgexc(true, true, 1, A, 1, B, 1, Q, 1, Z, 1, IFST, ILST, W, 0, INFO);
      chkxer('DTGEXC', infoc.INFOT, NOUT, LERR, OK);
      NT += 8;
    });

    test('DTGSEN', () {
      srnamc.SRNAMT = 'DTGSEN';
      infoc.INFOT = 1;
      dtgsen(-1, true, true, SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M,
          TOLA, TOLB, RCV, W, 1, IW, 1, INFO);
      chkxer('DTGSEN', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dtgsen(1, true, true, SEL, -1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M,
          TOLA, TOLB, RCV, W, 1, IW, 1, INFO);
      chkxer('DTGSEN', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dtgsen(1, true, true, SEL, 1, A, 0, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA,
          TOLB, RCV, W, 1, IW, 1, INFO);
      chkxer('DTGSEN', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dtgsen(1, true, true, SEL, 1, A, 1, B, 0, R1, R2, R3, Q, 1, Z, 1, M, TOLA,
          TOLB, RCV, W, 1, IW, 1, INFO);
      chkxer('DTGSEN', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 14;
      dtgsen(1, true, true, SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 0, Z, 1, M, TOLA,
          TOLB, RCV, W, 1, IW, 1, INFO);
      chkxer('DTGSEN', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 16;
      dtgsen(1, true, true, SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 0, M, TOLA,
          TOLB, RCV, W, 1, IW, 1, INFO);
      chkxer('DTGSEN', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 22;
      dtgsen(0, true, true, SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA,
          TOLB, RCV, W, 1, IW, 1, INFO);
      chkxer('DTGSEN', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 22;
      dtgsen(1, true, true, SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA,
          TOLB, RCV, W, 1, IW, 1, INFO);
      chkxer('DTGSEN', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 22;
      dtgsen(2, true, true, SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA,
          TOLB, RCV, W, 1, IW, 1, INFO);
      chkxer('DTGSEN', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 24;
      dtgsen(0, true, true, SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA,
          TOLB, RCV, W, 20, IW, 0, INFO);
      chkxer('DTGSEN', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 24;
      dtgsen(1, true, true, SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA,
          TOLB, RCV, W, 20, IW, 0, INFO);
      chkxer('DTGSEN', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 24;
      dtgsen(2, true, true, SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA,
          TOLB, RCV, W, 20, IW, 1, INFO);
      chkxer('DTGSEN', infoc.INFOT, NOUT, LERR, OK);
      NT += 12;
    });

    test('DTGSNA', () {
      srnamc.SRNAMT = 'DTGSNA';
      infoc.INFOT = 1;
      dtgsna('/', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW,
          INFO);
      chkxer('DTGSNA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dtgsna('B', '/', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW,
          INFO);
      chkxer('DTGSNA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dtgsna('B', 'A', SEL, -1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW,
          INFO);
      chkxer('DTGSNA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dtgsna('B', 'A', SEL, 1, A, 0, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW,
          INFO);
      chkxer('DTGSNA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dtgsna('B', 'A', SEL, 1, A, 1, B, 0, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW,
          INFO);
      chkxer('DTGSNA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dtgsna('E', 'A', SEL, 1, A, 1, B, 1, Q, 0, U, 1, R1, R2, 1, M, W, 1, IW,
          INFO);
      chkxer('DTGSNA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 12;
      dtgsna('E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 0, R1, R2, 1, M, W, 1, IW,
          INFO);
      chkxer('DTGSNA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 15;
      dtgsna('E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 0, M, W, 1, IW,
          INFO);
      chkxer('DTGSNA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 18;
      dtgsna('E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 0, IW,
          INFO);
      chkxer('DTGSNA', infoc.INFOT, NOUT, LERR, OK);
      NT += 9;
    });

    test('DTGSYL', () {
      srnamc.SRNAMT = 'DTGSYL';
      infoc.INFOT = 1;
      dtgsyl('/', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1,
          IW, INFO);
      chkxer('DTGSYL', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dtgsyl('N', -1, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W,
          1, IW, INFO);
      chkxer('DTGSYL', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dtgsyl('N', 0, 0, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1,
          IW, INFO);
      chkxer('DTGSYL', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dtgsyl('N', 0, 1, 0, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1,
          IW, INFO);
      chkxer('DTGSYL', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dtgsyl('N', 0, 1, 1, A, 0, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1,
          IW, INFO);
      chkxer('DTGSYL', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dtgsyl('N', 0, 1, 1, A, 1, B, 0, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1,
          IW, INFO);
      chkxer('DTGSYL', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dtgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 0, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1,
          IW, INFO);
      chkxer('DTGSYL', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 12;
      dtgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 0, V, 1, Z, 1, SCALE, DIF, W, 1,
          IW, INFO);
      chkxer('DTGSYL', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 14;
      dtgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 0, Z, 1, SCALE, DIF, W, 1,
          IW, INFO);
      chkxer('DTGSYL', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 16;
      dtgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 0, SCALE, DIF, W, 1,
          IW, INFO);
      chkxer('DTGSYL', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 20;
      dtgsyl('N', 1, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1,
          IW, INFO);
      chkxer('DTGSYL', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 20;
      dtgsyl('N', 2, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1,
          IW, INFO);
      chkxer('DTGSYL', infoc.INFOT, NOUT, LERR, OK);
      NT += 12;
    });
  }

  // Print a summary line.

  if (infoc.OK.value) {
    NOUT.println(
        ' ${PATH.a3} routines passed the tests of the error exits (${NT.i3} tests done)');
  } else {
    NOUT.println(
        ' *** ${PATH.a3} routines failed the tests of the error exits ***');
  }
}
